/**
 * QuPath Groovy Script: ROI 내 Point Annotation 분포 분석 (메모리 최적화 버전)
 * 
 * 대용량 데이터셋을 위한 최적화:
 * - KD-Tree를 사용한 효율적인 최근접 이웃 검색
 * - 대규모 데이터의 경우 샘플링 기반 분석
 * - 거리 행렬 대신 on-the-fly 계산
 * 
 * 분석 항목:
 * 1. 기본 통계 (개수, 밀도)
 * 2. 최근접 이웃 거리 (Nearest Neighbor Distance) 분석
 * 3. Clark-Evans Index (Clustering/Dispersion 지표)
 * 4. Ripley's K-function 기반 클러스터링 분석
 * 5. Quadrat 분석 (공간 분포 균일성)
 * 6. 세포 유형별 분석
 */

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.roi.PointsROI
import qupath.lib.roi.interfaces.ROI
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.index.kdtree.KdTree
import org.locationtech.jts.index.kdtree.KdNode

// ============================================================
// 설정
// ============================================================
def roiClassName = "ROI"  // ROI 클래스 이름
def pointClasses = ["lymphocyte", "connective-tissue-cell", "neutrophil", "plasma-cell", "eosinophil"]

// 분석 파라미터
def ripleyRadii = [10, 20, 30, 50, 75, 100, 150, 200] as double[]  // Ripley's K 분석 반경 (픽셀)
def quadratGridSize = 10  // Quadrat 분석용 그리드 크기 (n x n)
def neighborK = 5  // K번째 최근접 이웃까지 분석

// 샘플링 설정
def maxPointsForFullAnalysis = 10000  // 이 수 이하면 전체 분석
def sampleSizeForLargeData = 5000     // 대용량 데이터 샘플 크기
def ripleyMaxPoints = 3000            // Ripley 분석용 최대 포인트

// ============================================================
// KD-Tree 기반 유틸리티 클래스
// ============================================================
class SpatialAnalyzer {
    
    // 두 점 사이의 유클리드 거리 계산
    static double distance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2))
    }
    
    // KD-Tree 구축
    static KdTree buildKdTree(List<double[]> coords) {
        KdTree tree = new KdTree()
        coords.each { coord ->
            tree.insert(new Coordinate(coord[0], coord[1]))
        }
        return tree
    }
    
    // KD-Tree를 사용하여 K개의 최근접 이웃 찾기
    static List<Double> findKNearestDistances(KdTree tree, double x, double y, int k, List<double[]> allCoords) {
        // 초기 검색 반경 설정 (점 밀도 기반 추정)
        double searchRadius = 100  // 시작 반경
        List<Double> distances = []
        
        // 반경을 점점 늘려가며 검색
        while (distances.size() < k + 1 && searchRadius < 100000) {
            distances = []
            def envelope = new Envelope(x - searchRadius, x + searchRadius, y - searchRadius, y + searchRadius)
            def results = tree.query(envelope)
            
            results.each { node ->
                def coord = node.getCoordinate()
                double d = distance(x, y, coord.x, coord.y)
                if (d > 0) {  // 자기 자신 제외
                    distances << d
                }
            }
            searchRadius *= 2
        }
        
        distances.sort()
        return distances.take(k)
    }
    
    // 단일 최근접 이웃 거리 (효율적 버전)
    static double findNearestNeighborDistance(List<double[]> coords, int index, KdTree tree) {
        double x = coords[index][0]
        double y = coords[index][1]
        
        double searchRadius = 50
        double minDist = Double.MAX_VALUE
        
        while (minDist == Double.MAX_VALUE && searchRadius < 100000) {
            def envelope = new Envelope(x - searchRadius, x + searchRadius, y - searchRadius, y + searchRadius)
            def results = tree.query(envelope)
            
            results.each { node ->
                def coord = node.getCoordinate()
                double d = distance(x, y, coord.x, coord.y)
                if (d > 0 && d < minDist) {
                    minDist = d
                }
            }
            searchRadius *= 2
        }
        
        return minDist
    }
    
    // 모든 포인트의 최근접 이웃 거리 계산 (KD-Tree 사용)
    static double[] computeAllNNDistances(List<double[]> coords, KdTree tree) {
        double[] distances = new double[coords.size()]
        coords.eachWithIndex { coord, i ->
            distances[i] = findNearestNeighborDistance(coords, i, tree)
            if ((i + 1) % 10000 == 0) {
                println "  진행: ${i + 1}/${coords.size()} 포인트 처리됨"
            }
        }
        return distances
    }
    
    // 샘플링된 포인트의 K번째 최근접 이웃 거리
    static double[] computeSampledKthNNDistances(List<double[]> coords, KdTree tree, int k, int sampleSize) {
        def random = new Random(42)  // 재현성을 위한 시드
        def indices = (0..<coords.size()).toList()
        Collections.shuffle(indices, random)
        def sampleIndices = indices.take(Math.min(sampleSize, coords.size()))
        
        double[] distances = new double[sampleIndices.size()]
        sampleIndices.eachWithIndex { idx, i ->
            def coord = coords[idx]
            def kDistances = findKNearestDistances(tree, coord[0], coord[1], k, coords)
            distances[i] = kDistances.size() >= k ? kDistances[k - 1] : Double.NaN
            
            if ((i + 1) % 1000 == 0) {
                println "  진행: ${i + 1}/${sampleIndices.size()} 샘플 처리됨"
            }
        }
        return distances
    }
    
    // 기본 통계 계산
    static Map<String, Double> computeBasicStats(double[] values) {
        if (values.length == 0) {
            return [mean: Double.NaN, std: Double.NaN, min: Double.NaN, max: Double.NaN, median: Double.NaN]
        }
        
        def validValues = values.findAll { !Double.isNaN(it) && it != Double.MAX_VALUE }
        if (validValues.isEmpty()) {
            return [mean: Double.NaN, std: Double.NaN, min: Double.NaN, max: Double.NaN, median: Double.NaN]
        }
        
        double sum = validValues.sum()
        double mean = sum / validValues.size()
        double variance = validValues.collect { Math.pow(it - mean, 2) }.sum() / validValues.size()
        double std = Math.sqrt(variance)
        
        validValues.sort()
        double median = validValues.size() % 2 == 0 ? 
            (validValues[(int)(validValues.size()/2) - 1] + validValues[(int)(validValues.size()/2)]) / 2 :
            validValues[(int)(validValues.size()/2)]
        
        return [
            mean: mean,
            std: std,
            min: validValues.min(),
            max: validValues.max(),
            median: median,
            count: validValues.size()
        ]
    }
    
    // Clark-Evans Index 계산 (R)
    static double computeClarkEvansIndex(double meanNNDistance, double area, int n) {
        if (n < 2 || area <= 0) return Double.NaN
        double expectedMean = 0.5 * Math.sqrt(area / n)
        return meanNNDistance / expectedMean
    }
    
    // Clark-Evans Z-score 계산
    static double computeClarkEvansZScore(double meanNNDistance, double area, int n) {
        if (n < 2 || area <= 0) return Double.NaN
        double expectedMean = 0.5 * Math.sqrt(area / n)
        double se = 0.26136 / Math.sqrt(n * n / area)
        return (meanNNDistance - expectedMean) / se
    }
    
    // Ripley's K-function 계산 (샘플링 기반)
    static Map<Double, Double> computeRipleysK(List<double[]> coords, KdTree tree, double area, double[] radii, int maxPoints) {
        def random = new Random(42)
        def sampleCoords = coords
        
        // 포인트가 너무 많으면 샘플링
        if (coords.size() > maxPoints) {
            def indices = (0..<coords.size()).toList()
            Collections.shuffle(indices, random)
            sampleCoords = indices.take(maxPoints).collect { coords[it] }
            println "  Ripley's K: ${coords.size()}개 중 ${maxPoints}개 샘플링"
        }
        
        int n = sampleCoords.size()
        if (n < 2) return [:]
        
        // 전체 포인트에 대한 밀도 보정 계수
        double densityCorrection = coords.size() / (double) n
        
        Map<Double, Double> kValues = [:]
        
        radii.each { r ->
            println "  반경 ${r} 분석 중..."
            double sumK = 0
            
            sampleCoords.each { coord ->
                def envelope = new Envelope(coord[0] - r, coord[0] + r, coord[1] - r, coord[1] + r)
                def results = tree.query(envelope)
                
                results.each { node ->
                    def nodeCoord = node.getCoordinate()
                    double d = distance(coord[0], coord[1], nodeCoord.x, nodeCoord.y)
                    if (d > 0 && d <= r) {
                        sumK += 1
                    }
                }
            }
            
            // K(r) = A/n² * sum(I(d_ij <= r))
            // 샘플링 시 원래 포인트 수로 정규화
            double k = (area / (coords.size() * (coords.size() - 1))) * sumK * densityCorrection
            kValues[r] = k
        }
        
        return kValues
    }
    
    // Ripley's L-function 계산
    static Map<Double, Double> computeRipleysL(Map<Double, Double> kValues) {
        Map<Double, Double> lValues = [:]
        kValues.each { r, k ->
            double l = Math.sqrt(k / Math.PI) - r
            lValues[r] = l
        }
        return lValues
    }
    
    // Quadrat 분석
    static Map<String, Object> computeQuadratAnalysis(List<double[]> coords, double minX, double minY, 
                                                       double maxX, double maxY, int gridSize) {
        double width = maxX - minX
        double height = maxY - minY
        double cellWidth = width / gridSize
        double cellHeight = height / gridSize
        
        int[][] counts = new int[gridSize][gridSize]
        coords.each { coord ->
            int col = Math.min((int)((coord[0] - minX) / cellWidth), gridSize - 1)
            int row = Math.min((int)((coord[1] - minY) / cellHeight), gridSize - 1)
            if (col >= 0 && col < gridSize && row >= 0 && row < gridSize) {
                counts[row][col]++
            }
        }
        
        def allCounts = []
        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
                allCounts << counts[i][j]
            }
        }
        
        double mean = allCounts.sum() / allCounts.size()
        double variance = allCounts.collect { Math.pow(it - mean, 2) }.sum() / allCounts.size()
        double vmr = mean > 0 ? variance / mean : Double.NaN
        double chiSquare = mean > 0 ? allCounts.collect { Math.pow(it - mean, 2) / mean }.sum() : Double.NaN
        int df = allCounts.size() - 1
        
        return [
            counts: counts,
            mean: mean,
            variance: variance,
            vmr: vmr,
            chiSquare: chiSquare,
            degreesOfFreedom: df,
            interpretation: vmr < 0.9 ? "균일 분포" : (vmr > 1.1 ? "클러스터링" : "랜덤 분포")
        ]
    }
    
    // 히스토그램 생성
    static Map<String, Object> createHistogram(double[] values, int bins) {
        def validValues = values.findAll { !Double.isNaN(it) && it != Double.MAX_VALUE }
        if (validValues.isEmpty()) return [:]
        
        double min = validValues.min()
        double max = validValues.max()
        if (max == min) max = min + 1
        double binWidth = (max - min) / bins
        
        int[] histogram = new int[bins]
        validValues.each { v ->
            int bin = Math.min((int)((v - min) / binWidth), bins - 1)
            histogram[bin]++
        }
        
        def binEdges = (0..bins).collect { min + it * binWidth }
        
        return [
            histogram: histogram,
            binEdges: binEdges,
            binWidth: binWidth
        ]
    }
    
    // 랜덤 샘플링
    static List<double[]> sampleCoordinates(List<double[]> coords, int sampleSize, long seed = 42) {
        if (coords.size() <= sampleSize) return coords
        
        def random = new Random(seed)
        def indices = (0..<coords.size()).toList()
        Collections.shuffle(indices, random)
        return indices.take(sampleSize).collect { coords[it] }
    }
    
    // Cross-type 최근접 거리 (KD-Tree 사용)
    static double[] computeCrossTypeNNDistances(List<double[]> coordsA, List<double[]> coordsB, int maxSample = 1000) {
        if (coordsA.isEmpty() || coordsB.isEmpty()) return new double[0]
        
        // B의 KD-Tree 구축
        KdTree treeB = buildKdTree(coordsB)
        
        // A를 샘플링
        def sampleA = sampleCoordinates(coordsA, maxSample)
        double[] distances = new double[sampleA.size()]
        
        sampleA.eachWithIndex { coord, i ->
            double minDist = Double.MAX_VALUE
            double searchRadius = 100
            
            while (minDist == Double.MAX_VALUE && searchRadius < 100000) {
                def envelope = new Envelope(coord[0] - searchRadius, coord[0] + searchRadius, 
                                           coord[1] - searchRadius, coord[1] + searchRadius)
                def results = treeB.query(envelope)
                
                results.each { node ->
                    def nodeCoord = node.getCoordinate()
                    double d = distance(coord[0], coord[1], nodeCoord.x, nodeCoord.y)
                    if (d < minDist) minDist = d
                }
                searchRadius *= 2
            }
            distances[i] = minDist
        }
        
        return distances
    }
}

// ============================================================
// 메인 분석 코드
// ============================================================
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
def cal = server.getPixelCalibration()

// 픽셀 크기 (μm)
double pixelSizeUM = cal.getAveragedPixelSizeMicrons()
if (Double.isNaN(pixelSizeUM) || pixelSizeUM <= 0) {
    pixelSizeUM = 0.5
    println "경고: 픽셀 크기를 가져올 수 없어 기본값(${pixelSizeUM} μm)을 사용합니다."
}

println "=" * 80
println "ROI 내 Point Annotation 분포 분석 (메모리 최적화 버전)"
println "=" * 80
println "이미지: ${server.getMetadata().getName()}"
println "픽셀 크기: ${String.format('%.4f', pixelSizeUM)} μm"
println "=" * 80

// ROI 클래스 annotation 가져오기
def roiAnnotations = hierarchy.getAnnotationObjects().findAll { 
    it.getPathClass()?.getName() == roiClassName 
}

if (roiAnnotations.isEmpty()) {
    println "경고: '${roiClassName}' 클래스의 annotation을 찾을 수 없습니다."
    return
}

println "\n발견된 ROI annotation 수: ${roiAnnotations.size()}"

// 결과 저장용
def allResults = []

// 각 ROI에 대해 분석 수행
roiAnnotations.eachWithIndex { roiAnnotation, roiIndex ->
    println "\n" + "=" * 80
    println "ROI #${roiIndex + 1} 분석"
    println "=" * 80
    
    def roiGeometry = roiAnnotation.getROI().getGeometry()
    double roiAreaPixels = roiGeometry.getArea()
    double roiAreaUM2 = roiAreaPixels * pixelSizeUM * pixelSizeUM
    double roiAreaMM2 = roiAreaUM2 / 1e6
    
    def envelope = roiGeometry.getEnvelopeInternal()
    double minX = envelope.getMinX()
    double minY = envelope.getMinY()
    double maxX = envelope.getMaxX()
    double maxY = envelope.getMaxY()
    
    println "ROI 면적: ${String.format('%.2f', roiAreaUM2)} μm² (${String.format('%.4f', roiAreaMM2)} mm²)"
    println "ROI 경계: [${String.format('%.1f', minX)}, ${String.format('%.1f', minY)}] - [${String.format('%.1f', maxX)}, ${String.format('%.1f', maxY)}]"
    
    // ROI 내부의 Point annotation 찾기
    def gf = new GeometryFactory()
    def allPointAnnotations = hierarchy.getAnnotationObjects().findAll { ann ->
        def pathClass = ann.getPathClass()?.getName()
        if (pathClass == null || !pointClasses.contains(pathClass)) return false
        
        def annROI = ann.getROI()
        if (!(annROI instanceof PointsROI)) return false
        
        def points = annROI.getAllPoints()
        return points.any { pt ->
            def coord = new Coordinate(pt.getX(), pt.getY())
            def point = gf.createPoint(coord)
            roiGeometry.contains(point)
        }
    }
    
    if (allPointAnnotations.isEmpty()) {
        println "이 ROI 내에 Point annotation이 없습니다."
        return
    }
    
    // 전체 포인트 좌표 추출
    def allCoords = []
    def coordsByClass = [:]
    pointClasses.each { coordsByClass[it] = [] }
    
    allPointAnnotations.each { ann ->
        def className = ann.getPathClass()?.getName()
        def roi = ann.getROI()
        if (roi instanceof PointsROI) {
            roi.getAllPoints().each { pt ->
                def coord = new Coordinate(pt.getX(), pt.getY())
                def point = gf.createPoint(coord)
                if (roiGeometry.contains(point)) {
                    def coordArray = [pt.getX(), pt.getY()] as double[]
                    allCoords << coordArray
                    if (coordsByClass.containsKey(className)) {
                        coordsByClass[className] << coordArray
                    }
                }
            }
        }
    }
    
    int totalPoints = allCoords.size()
    println "\n총 포인트 수: ${totalPoints}"
    
    boolean useSampling = totalPoints > maxPointsForFullAnalysis
    if (useSampling) {
        println ">>> 대용량 데이터 감지: 샘플링 기반 분석 수행 (샘플 크기: ${sampleSizeForLargeData})"
    }
    
    // ============================================================
    // 1. 세포 유형별 기본 통계
    // ============================================================
    println "\n" + "-" * 40
    println "1. 세포 유형별 기본 통계"
    println "-" * 40
    
    def cellTypeCounts = [:]
    def cellTypeDensities = [:]
    
    pointClasses.each { className ->
        int count = coordsByClass[className].size()
        double density = count / roiAreaMM2
        cellTypeCounts[className] = count
        cellTypeDensities[className] = density
        
        if (count > 0) {
            println "${className}:"
            println "  - 개수: ${count}"
            println "  - 밀도: ${String.format('%.2f', density)} cells/mm²"
            println "  - 비율: ${String.format('%.1f', count * 100.0 / totalPoints)}%"
        }
    }
    
    println "\n전체 세포 밀도: ${String.format('%.2f', totalPoints / roiAreaMM2)} cells/mm²"
    
    // ============================================================
    // 2. KD-Tree 구축 및 최근접 이웃 거리 분석
    // ============================================================
    println "\n" + "-" * 40
    println "2. 최근접 이웃 거리 (Nearest Neighbor Distance) 분석"
    println "-" * 40
    
    if (totalPoints >= 2) {
        println "KD-Tree 구축 중..."
        def kdTree = SpatialAnalyzer.buildKdTree(allCoords)
        println "KD-Tree 구축 완료"
        
        // 분석용 좌표 (샘플링 또는 전체)
        def analysisCoords = useSampling ? 
            SpatialAnalyzer.sampleCoordinates(allCoords, sampleSizeForLargeData) : 
            allCoords
        
        if (useSampling) {
            println "샘플링: ${totalPoints}개 중 ${analysisCoords.size()}개 선택"
        }
        
        // 1차 최근접 이웃 거리 계산
        println "\n1차 최근접 이웃 거리 계산 중..."
        def nnDistances = new double[analysisCoords.size()]
        analysisCoords.eachWithIndex { coord, i ->
            double minDist = Double.MAX_VALUE
            double searchRadius = 50
            
            while (minDist == Double.MAX_VALUE && searchRadius < 100000) {
                def searchEnv = new Envelope(coord[0] - searchRadius, coord[0] + searchRadius, 
                                            coord[1] - searchRadius, coord[1] + searchRadius)
                def results = kdTree.query(searchEnv)
                
                results.each { node ->
                    def nodeCoord = node.getCoordinate()
                    double d = SpatialAnalyzer.distance(coord[0], coord[1], nodeCoord.x, nodeCoord.y)
                    if (d > 0 && d < minDist) minDist = d
                }
                searchRadius *= 2
            }
            nnDistances[i] = minDist
            
            if ((i + 1) % 5000 == 0) {
                println "  진행: ${i + 1}/${analysisCoords.size()}"
            }
        }
        
        def nnDistancesUM = nnDistances.collect { it * pixelSizeUM } as double[]
        def nnStats = SpatialAnalyzer.computeBasicStats(nnDistancesUM)
        
        println "\n1차 최근접 이웃 거리 (μm):"
        println "  - 평균: ${String.format('%.2f', nnStats.mean)}"
        println "  - 표준편차: ${String.format('%.2f', nnStats.std)}"
        println "  - 중앙값: ${String.format('%.2f', nnStats.median)}"
        println "  - 최소: ${String.format('%.2f', nnStats.min)}"
        println "  - 최대: ${String.format('%.2f', nnStats.max)}"
        
        // K번째 최근접 이웃 분석 (샘플링)
        if (neighborK > 1) {
            println "\n${neighborK}번째 최근접 이웃까지 분석 중 (샘플 기반)..."
            def kthSampleSize = Math.min(1000, analysisCoords.size())
            
            (2..neighborK).each { k ->
                def kthDistances = SpatialAnalyzer.computeSampledKthNNDistances(allCoords, kdTree, k, kthSampleSize)
                def kthDistancesUM = kthDistances.collect { it * pixelSizeUM } as double[]
                def kthStats = SpatialAnalyzer.computeBasicStats(kthDistancesUM)
                
                println "\n${k}번째 최근접 이웃 거리 (μm, n=${kthSampleSize} 샘플):"
                println "  - 평균: ${String.format('%.2f', kthStats.mean)}"
                println "  - 표준편차: ${String.format('%.2f', kthStats.std)}"
                println "  - 중앙값: ${String.format('%.2f', kthStats.median)}"
            }
        }
        
        // ============================================================
        // 3. Clark-Evans Index
        // ============================================================
        println "\n" + "-" * 40
        println "3. Clark-Evans Index (공간 분포 패턴)"
        println "-" * 40
        
        double ceIndex = SpatialAnalyzer.computeClarkEvansIndex(nnStats.mean, roiAreaUM2, totalPoints)
        double ceZScore = SpatialAnalyzer.computeClarkEvansZScore(nnStats.mean, roiAreaUM2, totalPoints)
        
        println "R index: ${String.format('%.4f', ceIndex)}"
        println "Z-score: ${String.format('%.4f', ceZScore)}"
        println "해석:"
        if (ceIndex < 0.9) {
            println "  → 클러스터링 경향 (R < 1)"
        } else if (ceIndex > 1.1) {
            println "  → 균일 분포 경향 (R > 1)"
        } else {
            println "  → 랜덤 분포에 가까움 (R ≈ 1)"
        }
        if (Math.abs(ceZScore) > 1.96) {
            println "  → 통계적으로 유의함 (|Z| > 1.96, p < 0.05)"
        } else {
            println "  → 통계적으로 유의하지 않음 (|Z| ≤ 1.96)"
        }
        
        // ============================================================
        // 4. Ripley's K/L-function 분석
        // ============================================================
        println "\n" + "-" * 40
        println "4. Ripley's K/L-function 분석 (다중 스케일 클러스터링)"
        println "-" * 40
        
        def kValues = SpatialAnalyzer.computeRipleysK(allCoords, kdTree, roiAreaPixels, ripleyRadii, ripleyMaxPoints)
        def lValues = SpatialAnalyzer.computeRipleysL(kValues)
        
        println "\n반경(μm)\tK(r)\t\tL(r)(μm)\t해석"
        println "-" * 60
        ripleyRadii.each { r ->
            double rUM = r * pixelSizeUM
            double k = kValues[r] ?: Double.NaN
            double l = lValues[r] ?: Double.NaN
            double lUM = l * pixelSizeUM
            String interpretation = lUM > 5 ? "클러스터링" : (lUM < -5 ? "균일" : "랜덤")
            println "${String.format('%.1f', rUM)}\t\t${String.format('%.2f', k)}\t\t${String.format('%.2f', lUM)}\t\t${interpretation}"
        }
        
        // ============================================================
        // 5. Quadrat 분석
        // ============================================================
        println "\n" + "-" * 40
        println "5. Quadrat 분석 (공간 균일성)"
        println "-" * 40
        
        def quadratResult = SpatialAnalyzer.computeQuadratAnalysis(allCoords, minX, minY, maxX, maxY, quadratGridSize)
        
        println "그리드 크기: ${quadratGridSize} x ${quadratGridSize}"
        println "Quadrat당 평균 포인트 수: ${String.format('%.2f', quadratResult.mean)}"
        println "분산: ${String.format('%.2f', quadratResult.variance)}"
        println "VMR (Variance-to-Mean Ratio): ${String.format('%.4f', quadratResult.vmr)}"
        println "Chi-square: ${String.format('%.2f', quadratResult.chiSquare)} (df=${quadratResult.degreesOfFreedom})"
        println "해석: ${quadratResult.interpretation}"
        
        // Quadrat 히트맵 (간소화 출력)
        println "\nQuadrat 카운트 분포 (min/max/mean):"
        def flatCounts = []
        quadratResult.counts.each { row -> flatCounts.addAll(row.toList()) }
        println "  최소: ${flatCounts.min()}, 최대: ${flatCounts.max()}, 평균: ${String.format('%.1f', flatCounts.sum() / flatCounts.size())}"
        
        // ============================================================
        // 6. 거리 분포 히스토그램
        // ============================================================
        println "\n" + "-" * 40
        println "6. 최근접 이웃 거리 분포 히스토그램"
        println "-" * 40
        
        def histResult = SpatialAnalyzer.createHistogram(nnDistancesUM, 10)
        if (!histResult.isEmpty()) {
            int maxCount = histResult.histogram.max()
            println "\n구간 (μm)\t\t빈도\t그래프"
            println "-" * 50
            histResult.histogram.eachWithIndex { count, i ->
                def low = histResult.binEdges[i]
                def high = histResult.binEdges[i + 1]
                int barLen = maxCount > 0 ? (int)(count * 30 / maxCount) : 0
                def bar = "█" * Math.max(1, barLen)
                println "${String.format('%6.1f', low)} - ${String.format('%6.1f', high)}\t${String.format('%6d', count)}\t${bar}"
            }
        }
        
        // ============================================================
        // 7. 세포 유형별 상세 분석
        // ============================================================
        println "\n" + "-" * 40
        println "7. 세포 유형별 상세 공간 분석"
        println "-" * 40
        
        pointClasses.each { className ->
            def classCoords = coordsByClass[className]
            if (classCoords.size() >= 10) {  // 최소 10개 이상
                println "\n>>> ${className} (n=${classCoords.size()}) <<<"
                
                def classTree = SpatialAnalyzer.buildKdTree(classCoords)
                def classSample = SpatialAnalyzer.sampleCoordinates(classCoords, Math.min(1000, classCoords.size()))
                
                double[] classNNDist = new double[classSample.size()]
                classSample.eachWithIndex { coord, i ->
                    double minDist = Double.MAX_VALUE
                    double searchRadius = 100
                    while (minDist == Double.MAX_VALUE && searchRadius < 100000) {
                        def searchEnv = new Envelope(coord[0] - searchRadius, coord[0] + searchRadius,
                                                    coord[1] - searchRadius, coord[1] + searchRadius)
                        def results = classTree.query(searchEnv)
                        results.each { node ->
                            def nodeCoord = node.getCoordinate()
                            double d = SpatialAnalyzer.distance(coord[0], coord[1], nodeCoord.x, nodeCoord.y)
                            if (d > 0 && d < minDist) minDist = d
                        }
                        searchRadius *= 2
                    }
                    classNNDist[i] = minDist * pixelSizeUM
                }
                
                def classStats = SpatialAnalyzer.computeBasicStats(classNNDist)
                println "1차 최근접 이웃 거리 (μm): 평균=${String.format('%.2f', classStats.mean)}, SD=${String.format('%.2f', classStats.std)}"
                
                double classCE = SpatialAnalyzer.computeClarkEvansIndex(classStats.mean, roiAreaUM2, classCoords.size())
                println "Clark-Evans Index: ${String.format('%.4f', classCE)}"
                
                if (classCE < 0.9) {
                    println "  → 클러스터링 패턴"
                } else if (classCE > 1.1) {
                    println "  → 균일 분포 패턴"
                } else {
                    println "  → 랜덤 분포"
                }
            } else if (classCoords.size() > 0) {
                println "\n>>> ${className}: n=${classCoords.size()} (분석하기에 불충분, 최소 10개 필요)"
            }
        }
        
        // ============================================================
        // 8. 세포 유형 간 상호작용 분석
        // ============================================================
        println "\n" + "-" * 40
        println "8. 세포 유형 간 최근접 거리 분석 (Cross-type Analysis)"
        println "-" * 40
        
        def activeClasses = pointClasses.findAll { coordsByClass[it].size() >= 10 }
        
        if (activeClasses.size() >= 2) {
            println "\n[유형A → 유형B] = A의 각 세포에서 가장 가까운 B 세포까지의 평균 거리"
            println ""
            
            activeClasses.each { classA ->
                activeClasses.each { classB ->
                    if (classA != classB) {
                        def crossDistances = SpatialAnalyzer.computeCrossTypeNNDistances(
                            coordsByClass[classA], coordsByClass[classB], 500)
                        
                        def crossDistancesUM = crossDistances.collect { it * pixelSizeUM } as double[]
                        def crossStats = SpatialAnalyzer.computeBasicStats(crossDistancesUM)
                        
                        println "[${classA} → ${classB}]: ${String.format('%.2f', crossStats.mean)} ± ${String.format('%.2f', crossStats.std)} μm"
                    }
                }
            }
        } else {
            println "충분한 세포 유형이 없어 상호작용 분석을 수행할 수 없습니다. (최소 2개 유형, 각 10개 이상)"
        }
        
        // 결과 저장
        allResults << [
            roiIndex: roiIndex + 1,
            roiArea_mm2: roiAreaMM2,
            totalPoints: totalPoints,
            density_cells_per_mm2: totalPoints / roiAreaMM2,
            cellTypeCounts: cellTypeCounts,
            cellTypeDensities: cellTypeDensities,
            nn1_mean_um: nnStats.mean,
            nn1_std_um: nnStats.std,
            clarkEvansIndex: ceIndex,
            clarkEvansZScore: ceZScore,
            quadratVMR: quadratResult.vmr
        ]
        
    } else {
        println "포인트 수가 부족하여 거리 분석을 수행할 수 없습니다. (최소 2개 필요)"
    }
}

// ============================================================
// 요약 출력
// ============================================================
println "\n\n" + "=" * 80
println "분석 요약"
println "=" * 80

if (allResults.size() > 0) {
    println "\nROI\t면적(mm²)\t총세포\t밀도(/mm²)\tNN1(μm)\t\tClark-Evans\tVMR\t분포"
    println "-" * 110
    
    allResults.each { r ->
        String pattern = ""
        if (r.clarkEvansIndex < 0.9) pattern = "클러스터링"
        else if (r.clarkEvansIndex > 1.1) pattern = "균일"
        else pattern = "랜덤"
        
        println "#${r.roiIndex}\t${String.format('%.4f', r.roiArea_mm2)}\t\t${r.totalPoints}\t${String.format('%.1f', r.density_cells_per_mm2)}\t\t${String.format('%.2f±%.2f', r.nn1_mean_um, r.nn1_std_um)}\t${String.format('%.3f', r.clarkEvansIndex)}\t\t${String.format('%.2f', r.quadratVMR)}\t${pattern}"
    }
    
    println "\n세포 유형별 총계:"
    pointClasses.each { className ->
        int total = allResults.collect { it.cellTypeCounts[className] ?: 0 }.sum()
        double avgDensity = allResults.collect { it.cellTypeDensities[className] ?: 0 }.sum() / allResults.size()
        if (total > 0) {
            println "  ${className}: ${total} (평균 밀도: ${String.format('%.2f', avgDensity)} cells/mm²)"
        }
    }
}

println "\n분석 완료!"
println "=" * 80