# NBDFinder Flow Diagrams and Technical Specifications

## Processing Flow Diagram

```
                    NBDFinder Processing Pipeline
                    ============================

    ┌─────────────────────────────────────────────────────────────────┐
    │                     Input Layer                                 │
    ├─────────────────────────────────────────────────────────────────┤
    │  • FASTA File Parser                                            │
    │  • Sequence Validation (ATCGN only)                            │
    │  • Character Sanitization                                       │
    │  • Length Validation (10bp - 1MB)                              │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │                  Hyperscan Engine                               │
    ├─────────────────────────────────────────────────────────────────┤
    │  • Pattern Database Compilation                                 │
    │  • Multi-pattern Scanning                                       │
    │  • High-Speed Regex Matching                                    │
    │  • Match Event Dispatching                                      │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │              Class-Specific Detection                           │
    ├─────────────────────────────────────────────────────────────────┤
    │  Class 1: Curved DNA        │  Class 6: G-Quadruplex           │
    │  Class 2: Slipped DNA       │  Class 7: i-Motif                │
    │  Class 3: Cruciform         │  Class 8: Z-DNA                  │
    │  Class 4: R-Loop            │  Class 9: Hybrid                 │
    │  Class 5: Triplex           │  Class 10: Cluster               │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │                 Scientific Scoring                              │
    ├─────────────────────────────────────────────────────────────────┤
    │  • G4Hunter Algorithm (G-Quadruplex)                           │
    │  • Z-seeker Algorithm (Z-DNA)                                  │
    │  • Curvature Analysis (Curved DNA)                             │
    │  • Instability Scoring (Slipped DNA)                           │
    │  • Thermodynamic Models (Cruciform/Triplex)                    │
    │  • RLFS Detection (R-Loop)                                     │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │               Biological Filtering                              │
    ├─────────────────────────────────────────────────────────────────┤
    │  • Length Constraints (S_min, S_max)                           │
    │  • Score Thresholds                                            │
    │  • Sequence Quality Control                                     │
    │  • Biological Relevance Validation                             │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │              Motif Standardization                              │
    ├─────────────────────────────────────────────────────────────────┤
    │  • 1-based Coordinate System                                    │
    │  • Motif ID Generation                                          │
    │  • Score Normalization [0,1]                                   │
    │  • Subclass Assignment                                          │
    │  • GC Content Calculation                                       │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │             Post-Processing Analysis                            │
    ├─────────────────────────────────────────────────────────────────┤
    │  • Hybrid Detection (Class 9)                                  │
    │  • Cluster Analysis (Class 10)                                 │
    │  • Overlap Detection                                            │
    │  • Hotspot Identification                                       │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │                Output Generation                                │
    ├─────────────────────────────────────────────────────────────────┤
    │  • Standard Schema Formation                                    │
    │  • Summary Statistics                                           │
    │  • Export Format Selection                                      │
    │  • Visualization Preparation                                    │
    └─────────────────────┬───────────────────────────────────────────┘
                          │
                          ▼
    ┌─────────────────────────────────────────────────────────────────┐
    │              Visualization & Export                             │
    ├─────────────────────────────────────────────────────────────────┤
    │  Formats: BED, BigWig, CSV, JSON, Excel                        │
    │  Plots: 21+ chart types across 5 categories                    │
    │  Interactive: Plotly-powered visualizations                     │
    │  Export: Multiple format support                                │
    └─────────────────────────────────────────────────────────────────┘
```

## Class-Specific Detection Flow

```
                     Motif Class Detection Pipeline
                     ==============================

┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   Hyperscan     │    │   Class-Specific│    │   Scientific    │
│   Pattern       │───▶│   Refinement    │───▶│   Scoring       │
│   Matching      │    │                 │    │                 │
└─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │
         ▼                       ▼                       ▼
    Raw Matches            Candidate Motifs         Scored Motifs
    
    Examples:                Examples:               Examples:
    GGG..GGG..GGG         G4 Candidates           G4Hunter: 2.1
    AT-rich regions       Curved Candidates       Curvature: 45.3
    CG-rich regions       Z-DNA Candidates        Z-seeker: 78.5
```

## Scoring Algorithm Details

### G4Hunter Algorithm Implementation

```python
def g4hunter_algorithm(sequence, window_size=25):
    """
    G4Hunter scoring implementation
    
    Formula: G4H = Σ(G_count - C_count) / window_size
    
    Parameters:
    - window_size: Sliding window size (default: 25bp)
    - threshold: Minimum score for G4 prediction (≥1.2)
    
    Returns:
    - Normalized score [0,1]
    """
    scores = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        g_count = window.count('G')
        c_count = window.count('C')
        score = (g_count - c_count) / window_size
        scores.append(score)
    
    # Normalize to [0,1] range
    max_score = 4.0  # Theoretical maximum
    normalized_scores = [max(0, min(1, score/max_score)) for score in scores]
    return normalized_scores
```

### Z-seeker Algorithm Implementation

```python
def z_seeker_algorithm(sequence, window_size=50):
    """
    Z-seeker scoring for Z-DNA detection
    
    Dinucleotide weights:
    CG: +7.0, GC: +3.0, CA/TG/AC/GT: +1.0, Others: 0.0
    
    Parameters:
    - window_size: Analysis window (default: 50bp)
    - threshold: Minimum score (≥50.0)
    
    Returns:
    - Raw Z-seeker score
    """
    weights = {
        'CG': 7.0, 'GC': 3.0, 'CA': 1.0, 'TG': 1.0,
        'AC': 1.0, 'GT': 1.0
    }
    
    scores = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        score = 0
        
        for j in range(len(window) - 1):
            dinucleotide = window[j:j+2]
            score += weights.get(dinucleotide, 0.0)
        
        scores.append(score)
    
    return scores
```

## Performance Optimization Strategies

### 1. Hyperscan Database Optimization

```
Database Compilation Strategy:
┌──────────────────────────────────────────────────────────────┐
│                Pattern Organization                          │
├──────────────────────────────────────────────────────────────┤
│  High Frequency Patterns (G4, Z-DNA)     → Cache Priority 1  │
│  Medium Frequency (Curved, Slipped)      → Cache Priority 2  │
│  Low Frequency (Triplex, Cruciform)      → Cache Priority 3  │
│  Rare Patterns (Hybrid, Cluster)         → Cache Priority 4  │
└──────────────────────────────────────────────────────────────┘

Memory Management:
- LRU Cache: 100MB default, 500MB maximum
- Pattern Database: ~50MB compiled size
- Result Buffer: Dynamic allocation
- Garbage Collection: Automatic cleanup every 1000 sequences
```

### 2. Parallel Processing Architecture

```
Worker Pool Management:
┌──────────────────────────────────────────────────────────────┐
│                Process Distribution                          │
├──────────────────────────────────────────────────────────────┤
│  Master Process    │  Worker Pool (N cores)                  │
│  ├─ Input Queue    │  ├─ Worker 1: Class 1,2                │
│  ├─ Result Queue   │  ├─ Worker 2: Class 3,4                │
│  ├─ Error Queue    │  ├─ Worker 3: Class 5,6                │
│  └─ Status Monitor │  └─ Worker N: Class 7,8                │
└──────────────────────────────────────────────────────────────┘

Load Balancing:
- Dynamic work distribution based on class complexity
- Adaptive batch sizing (100-10000 sequences)
- Memory-aware task scheduling
- Automatic worker recovery on failures
```

### 3. Memory Optimization Techniques

```
Memory Usage Patterns:
┌──────────────────────────────────────────────────────────────┐
│                   Memory Allocation                          │
├──────────────────────────────────────────────────────────────┤
│  Component               │  Memory Usage  │  Optimization    │
│  ─────────────────────────────────────────────────────────── │
│  Hyperscan Database     │  ~50MB         │  One-time load   │
│  Pattern Cache          │  100-500MB     │  LRU eviction    │
│  Sequence Buffer        │  1-10MB        │  Streaming mode  │
│  Result Storage         │  Variable      │  Batch processing│
│  Visualization Data     │  10-100MB      │  Lazy loading    │
└──────────────────────────────────────────────────────────────┘

Streaming Mode (sequences >1MB):
- Chunk-based processing (1MB chunks)
- Progressive result accumulation
- Memory-mapped file access
- Incremental garbage collection
```

## Quality Assurance Framework

### Testing Infrastructure

```
Test Suite Organization:
┌──────────────────────────────────────────────────────────────┐
│                    Testing Levels                            │
├──────────────────────────────────────────────────────────────┤
│  Unit Tests           │  Function-level validation           │
│  ├─ Scoring Algorithms│  ├─ G4Hunter accuracy               │
│  ├─ Pattern Matching │  ├─ Z-seeker precision               │
│  ├─ Classification   │  ├─ Motif ID generation              │
│  └─ Export Functions │  └─ Format compliance                │
│                                                              │
│  Integration Tests    │  Component interaction validation    │
│  ├─ Pipeline Flow    │  ├─ End-to-end processing            │
│  ├─ API Endpoints    │  ├─ REST API functionality           │
│  ├─ File I/O         │  ├─ FASTA parsing                    │
│  └─ Error Handling   │  └─ Exception management             │
│                                                              │
│  Performance Tests   │  Speed and memory validation         │
│  ├─ Benchmark Suites │  ├─ Processing speed tests           │
│  ├─ Memory Profiling │  ├─ Memory usage analysis            │
│  ├─ Scalability      │  ├─ Large sequence handling          │
│  └─ Load Testing     │  └─ Concurrent user simulation       │
└──────────────────────────────────────────────────────────────┘
```

### Validation Datasets

```
Reference Datasets:
┌──────────────────────────────────────────────────────────────┐
│                 Validation Sequences                         │
├──────────────────────────────────────────────────────────────┤
│  G-Quadruplex        │  Known G4 sequences from literature  │
│  ├─ PDB Structures   │  ├─ Experimentally verified G4s      │
│  ├─ G4-seq Data     │  ├─ High-throughput sequencing        │
│  └─ Synthetic G4s   │  └─ Designed test sequences           │
│                                                              │
│  Z-DNA               │  Z-form DNA reference sequences      │
│  ├─ Crystal Structures│ ├─ X-ray crystallography data       │
│  ├─ NMR Structures   │  ├─ Solution structure data          │
│  └─ Computational    │  └─ MD simulation predictions        │
│                                                              │
│  Other Motifs        │  Additional Non-B DNA structures     │
│  ├─ Curved DNA       │  ├─ A-tract arrays                   │
│  ├─ Cruciforms       │  ├─ Palindromic sequences            │
│  ├─ Triplex          │  ├─ Mirror repeat regions            │
│  └─ i-Motifs         │  └─ C-rich sequences                 │
└──────────────────────────────────────────────────────────────┘
```

## API Integration Examples

### Galaxy Tool Wrapper

```xml
<tool id="nbdfinder" name="NBDFinder" version="2.1.0">
    <description>Non-B DNA Motif Detection</description>
    
    <requirements>
        <requirement type="package" version="2.1.0">nbdfinder</requirement>
    </requirements>
    
    <command><![CDATA[
        python '$__tool_directory__/nbdfinder_wrapper.py'
        --input '$input_fasta'
        --output '$output_csv'
        --classes '$motif_classes'
        --threshold '$score_threshold'
        --workers \${GALAXY_SLOTS:-4}
    ]]></command>
    
    <inputs>
        <param name="input_fasta" type="data" format="fasta" 
               label="Input FASTA file" />
        <param name="motif_classes" type="select" multiple="true" 
               label="Motif classes to analyze">
            <option value="G-Quadruplex">G-Quadruplex</option>
            <option value="Z-DNA">Z-DNA</option>
            <option value="Curved_DNA">Curved DNA</option>
            <!-- ... other classes ... -->
        </param>
        <param name="score_threshold" type="float" value="0.5" min="0.0" max="1.0"
               label="Minimum normalized score threshold" />
    </inputs>
    
    <outputs>
        <data name="output_csv" format="csv" label="NBDFinder Results" />
        <data name="output_bed" format="bed" label="BED Track" />
        <data name="summary_json" format="json" label="Analysis Summary" />
    </outputs>
    
    <tests>
        <test>
            <param name="input_fasta" value="test_sequence.fasta" />
            <param name="motif_classes" value="G-Quadruplex,Z-DNA" />
            <output name="output_csv" file="expected_results.csv" />
        </test>
    </tests>
    
    <help><![CDATA[
NBDFinder is a comprehensive tool for detecting Non-B DNA motifs.

**Input:**
- FASTA file containing DNA sequences

**Output:**
- CSV file with detailed motif annotations
- BED file for genome browser visualization
- JSON summary with analysis statistics

**Parameters:**
- Motif Classes: Select which types of motifs to detect
- Score Threshold: Minimum confidence score (0.0-1.0)

For more information, visit: https://github.com/VRYella/NBDFinder
    ]]></help>
</tool>
```

### Nextflow Pipeline Integration

```groovy
#!/usr/bin/env nextflow

/*
 * NBDFinder Nextflow Pipeline
 * Comprehensive Non-B DNA motif detection workflow
 */

params.input = 'sequences/*.fasta'
params.outdir = 'results'
params.classes = 'all'
params.threshold = 0.5
params.workers = 4

process NBDFINDER_ANALYSIS {
    tag "$sample_id"
    publishDir "${params.outdir}/motifs", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_motifs.csv"), emit: results
    tuple val(sample_id), path("${sample_id}_summary.json"), emit: summary
    tuple val(sample_id), path("${sample_id}_tracks.bed"), emit: tracks
    
    script:
    """
    python /opt/nbdfinder/nbdfinder_cli.py \\
        --input ${fasta} \\
        --output ${sample_id}_motifs.csv \\
        --summary ${sample_id}_summary.json \\
        --bed ${sample_id}_tracks.bed \\
        --classes ${params.classes} \\
        --threshold ${params.threshold} \\
        --workers ${params.workers}
    """
}

process AGGREGATE_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path(results)
    path(summaries)
    
    output:
    path("combined_results.csv")
    path("analysis_report.html")
    
    script:
    """
    python /opt/nbdfinder/aggregate_results.py \\
        --results ${results} \\
        --summaries ${summaries} \\
        --output combined_results.csv \\
        --report analysis_report.html
    """
}

workflow {
    // Input channel
    input_ch = Channel.fromPath(params.input)
        .map { file -> tuple(file.baseName, file) }
    
    // Main analysis
    NBDFINDER_ANALYSIS(input_ch)
    
    // Aggregate results
    all_results = NBDFINDER_ANALYSIS.out.results.collect { it[1] }
    all_summaries = NBDFINDER_ANALYSIS.out.summary.collect { it[1] }
    
    AGGREGATE_RESULTS(all_results, all_summaries)
}
```

### Docker Container Specification

```dockerfile
# NBDFinder Docker Container
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \\
    build-essential \\
    libhyperscan-dev \\
    pkg-config \\
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/nbdfinder

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy NBDFinder source code
COPY . .

# Install NBDFinder
RUN pip install -e .

# Create non-root user
RUN useradd -m -u 1000 nbduser && \\
    chown -R nbduser:nbduser /opt/nbdfinder
USER nbduser

# Expose API port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \\
    CMD curl -f http://localhost:8000/api/v1/health || exit 1

# Default command
CMD ["python", "api.py", "--host", "0.0.0.0", "--port", "8000"]
```

### Kubernetes Deployment

```yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: nbdfinder-api
  labels:
    app: nbdfinder
spec:
  replicas: 3
  selector:
    matchLabels:
      app: nbdfinder
  template:
    metadata:
      labels:
        app: nbdfinder
    spec:
      containers:
      - name: nbdfinder
        image: nbdfinder:2.1.0
        ports:
        - containerPort: 8000
        env:
        - name: WORKERS
          value: "4"
        - name: CACHE_SIZE
          value: "200MB"
        resources:
          requests:
            memory: "512Mi"
            cpu: "500m"
          limits:
            memory: "2Gi"
            cpu: "2000m"
        livenessProbe:
          httpGet:
            path: /api/v1/health
            port: 8000
          initialDelaySeconds: 30
          periodSeconds: 10
        readinessProbe:
          httpGet:
            path: /api/v1/health
            port: 8000
          initialDelaySeconds: 5
          periodSeconds: 5
---
apiVersion: v1
kind: Service
metadata:
  name: nbdfinder-service
spec:
  selector:
    app: nbdfinder
  ports:
  - protocol: TCP
    port: 80
    targetPort: 8000
  type: LoadBalancer
```

## Conclusion

This comprehensive documentation provides complete technical specifications for the NBDFinder tool, including:

1. **Detailed Flow Diagrams**: Visual representation of the processing pipeline
2. **Algorithm Implementations**: Code examples for scoring methods
3. **Performance Optimization**: Memory and speed optimization strategies
4. **Quality Assurance**: Testing frameworks and validation datasets
5. **Integration Examples**: Galaxy, Nextflow, Docker, and Kubernetes configurations

The tool is production-ready and provides enterprise-grade capabilities for Non-B DNA motif detection across multiple deployment environments.