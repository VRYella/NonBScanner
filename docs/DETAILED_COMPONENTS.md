# 🧬 NonBScanner Detailed Component Diagrams

## 📋 Table of Contents
1. [Motif Detection Classes Hierarchy](#motif-detection-classes-hierarchy)
2. [Hyperscan Pattern Matching Flow](#hyperscan-pattern-matching-flow)
3. [Visualization Pipeline](#visualization-pipeline)
4. [Function & Library Dependencies](#function--library-dependencies)
5. [Scoring Algorithm Workflow](#scoring-algorithm-workflow)
6. [Configuration & Class Management](#configuration--class-management)

## 🗂️ Motif Detection Classes Hierarchy

```mermaid
mindmap
  root((NonBScanner<br/>Motif Classes))
    Class 1: Curved DNA
      Global Curvature
      Local Curvature
      A-tract Phasing
    Class 2: Slipped DNA
      Direct Repeat
      STR (Short Tandem Repeat)
      Microsatellites
    Class 3: Cruciform
      Palindromic Inverted Repeat
      Four-way Junction
      Stem-Loop Formation
    Class 4: R-Loop
      R-Loop Formation Site
      RNA-DNA Hybrid
      Transcription Bubble
    Class 5: Triplex
      Mirror Repeat
      Sticky DNA
      Triple Helix
    Class 6: G-Quadruplex
      Canonical G4
      Bulged G4
      Relaxed G4
      Bipartite G4
      Multimeric G4
      Imperfect G4
      G-Triplex
    Class 7: i-Motif
      Canonical i-Motif
      Extended i-Motif
      AC-Motif
    Class 8: Z-DNA
      Classic Z-DNA
      eGZ (Extended GZ)
      Left-handed Helix
    Class 9: A-philic DNA
      A-tract Affinity
      AT-rich Regions
      Protein Binding
    Class 10: Hybrid
      Multi-class Overlap
      Composite Structures
      Complex Formations
    Class 11: Clusters
      Motif Hotspot
      Mixed Cluster
      High-density Region
```

## ⚡ Hyperscan Pattern Matching Flow

```mermaid
flowchart LR
    subgraph "Pattern Registry"
        REGEX_REG[🔍 Regex Registry<br/>core/regex_registry.py]
        PATTERNS[📝 Pattern Definitions<br/>Class-specific regex]
        COMPILE[⚙️ Pattern Compilation<br/>Hyperscan format]
    end
    
    subgraph "Hyperscan Engine"
        HS_DB[🚀 Hyperscan Database<br/>Compiled patterns]
        HS_SCAN[⚡ Fast Scanning<br/>Intel SSE/AVX]
        HS_MATCH[🎯 Match Callback<br/>Position tracking]
    end
    
    subgraph "Post-Processing"
        VALIDATE_MATCH[✅ Match Validation<br/>Length, context checks]
        SCORE_CALC[📊 Score Calculation<br/>Class-specific algorithms]
        QUALITY_FILTER[🎛️ Quality Filter<br/>Threshold application]
    end
    
    REGEX_REG --> PATTERNS
    PATTERNS --> COMPILE
    COMPILE --> HS_DB
    
    HS_DB --> HS_SCAN
    HS_SCAN --> HS_MATCH
    
    HS_MATCH --> VALIDATE_MATCH
    VALIDATE_MATCH --> SCORE_CALC
    SCORE_CALC --> QUALITY_FILTER
    
    QUALITY_FILTER --> RESULTS[📋 Motif Results]
```

## 📈 Visualization Pipeline

```mermaid
graph TD
    subgraph "Input Data"
        MOTIF_DATA[🧬 Motif Detection Results]
        SEQ_DATA[📄 Sequence Information]
        STATS_DATA[📊 Statistical Data]
    end
    
    subgraph "Plotting Engines"
        MATPLOTLIB[📊 Matplotlib<br/>Static plots]
        PLOTLY[📈 Plotly<br/>Interactive plots]
        SEABORN[🎨 Seaborn<br/>Statistical plots]
    end
    
    subgraph "Chart Types"
        DIST[📊 Distribution Plots<br/>Motif counts, lengths]
        HEAT[🔥 Heatmaps<br/>Position matrices]
        TRACK[📍 Track Plots<br/>Genome browser style]
        SUNBURST[☀️ Sunburst<br/>Hierarchical data]
        NETWORK[🕸️ Network<br/>Motif interactions]
        SCATTER[📈 Scatter<br/>Score correlations]
        BOX[📦 Box Plots<br/>Score distributions]
        VIOLIN[🎻 Violin Plots<br/>Density curves]
        MANHATTAN[🏙️ Manhattan<br/>Position-score plot]
        TSNE[🧬 t-SNE<br/>Dimensionality reduction]
        TREEMAP[🌳 Treemap<br/>Hierarchical proportions]
    end
    
    subgraph "Export Formats"
        PNG[🖼️ PNG/JPEG<br/>Static images]
        SVG[📐 SVG<br/>Vector graphics]
        HTML[🌐 HTML<br/>Interactive plots]
        PDF[📄 PDF<br/>Publication quality]
    end
    
    subgraph "Browser Integration"
        UCSC[🌍 UCSC Hub<br/>Track hub creation]
        IGV[🔬 IGV Session<br/>Session files]
        JBROWSE[🧬 JBrowse<br/>Config generation]
        BIGWIG[📊 BigWig<br/>Signal tracks]
        BED[📋 BED<br/>Feature annotation]
    end
    
    MOTIF_DATA --> MATPLOTLIB
    SEQ_DATA --> MATPLOTLIB
    STATS_DATA --> MATPLOTLIB
    
    MOTIF_DATA --> PLOTLY
    SEQ_DATA --> PLOTLY
    STATS_DATA --> PLOTLY
    
    MOTIF_DATA --> SEABORN
    STATS_DATA --> SEABORN
    
    MATPLOTLIB --> DIST
    MATPLOTLIB --> HEAT
    MATPLOTLIB --> BOX
    MATPLOTLIB --> VIOLIN
    
    PLOTLY --> TRACK
    PLOTLY --> SUNBURST
    PLOTLY --> NETWORK
    PLOTLY --> SCATTER
    PLOTLY --> MANHATTAN
    PLOTLY --> TREEMAP
    
    SEABORN --> TSNE
    
    DIST --> PNG
    HEAT --> PNG
    TRACK --> HTML
    SUNBURST --> HTML
    
    PNG --> PDF
    HTML --> SVG
    
    MOTIF_DATA --> UCSC
    MOTIF_DATA --> IGV
    MOTIF_DATA --> JBROWSE
    MOTIF_DATA --> BIGWIG
    MOTIF_DATA --> BED
```

## 🔧 Function & Library Dependencies

```mermaid
graph TB
    subgraph "Core Libraries"
        NUMPY[🔢 NumPy<br/>Numerical computing]
        PANDAS[🐼 Pandas<br/>Data manipulation]
        SCIPY[🧮 SciPy<br/>Scientific computing]
        BIOPYTHON[🧬 BioPython<br/>Sequence analysis]
        HYPERSCAN[🚀 Hyperscan<br/>Pattern matching]
    end
    
    subgraph "Web Frameworks"
        STREAMLIT[🌐 Streamlit<br/>Web interface]
        FASTAPI[⚡ FastAPI<br/>REST API]
        UVICORN[🦄 Uvicorn<br/>ASGI server]
        PYDANTIC[📝 Pydantic<br/>Data validation]
    end
    
    subgraph "Visualization"
        MATPLOTLIB[📊 Matplotlib<br/>Plotting]
        PLOTLY[📈 Plotly<br/>Interactive viz]
        SEABORN[🎨 Seaborn<br/>Statistical plots]
        VENN[⭕ Matplotlib-venn<br/>Venn diagrams]
        NETWORKX[🕸️ NetworkX<br/>Graph analysis]
    end
    
    subgraph "Processing & Performance"
        NUMBA[⚡ Numba<br/>JIT compilation]
        MULTIPROCESSING[🔄 Multiprocessing<br/>Parallel execution]
        CONCURRENT[🔀 Concurrent.futures<br/>Thread/process pools]
        FUNCTOOLS[🛠️ Functools<br/>Function utilities]
    end
    
    subgraph "File I/O & Export"
        OPENPYXL[📊 OpenPyXL<br/>Excel files]
        XLSXWRITER[✍️ XlsxWriter<br/>Excel writing]
        JSON[📄 JSON<br/>Data interchange]
        CSV[📋 CSV<br/>Tabular data]
        PATHLIB[📁 Pathlib<br/>File paths]
    end
    
    subgraph "Utility & Analysis"
        RE[🔍 Regex<br/>Pattern matching]
        COLLECTIONS[📦 Collections<br/>Data structures]
        ITERTOOLS[🔄 Itertools<br/>Iterator utilities]
        STATISTICS[📈 Statistics<br/>Statistical functions]
        RANDOM[🎲 Random<br/>Random generation]
        DATETIME[📅 Datetime<br/>Time handling]
    end
    
    subgraph "Scientific Computing"
        SKLEARN[🤖 Scikit-learn<br/>Machine learning]
        TSNE[🧬 t-SNE<br/>Dimensionality reduction]
        CLUSTERING[🎯 Clustering<br/>Data grouping]
    end
    
    NUMPY --> PANDAS
    NUMPY --> SCIPY
    NUMPY --> NUMBA
    
    PANDAS --> MATPLOTLIB
    PANDAS --> PLOTLY
    PANDAS --> SEABORN
    
    BIOPYTHON --> PANDAS
    
    STREAMLIT --> PANDAS
    STREAMLIT --> MATPLOTLIB
    STREAMLIT --> PLOTLY
    
    FASTAPI --> PYDANTIC
    FASTAPI --> UVICORN
    
    SKLEARN --> NUMPY
    SKLEARN --> PANDAS
    
    HYPERSCAN --> RE
    MULTIPROCESSING --> CONCURRENT
```

## ⚖️ Scoring Algorithm Workflow

```mermaid
stateDiagram-v2
    [*] --> SequenceInput
    
    SequenceInput --> PatternMatch : Hyperscan Detection
    
    PatternMatch --> G4Hunter : G-Quadruplex
    PatternMatch --> ZSeeker : Z-DNA
    PatternMatch --> CurvatureAnalysis : Curved DNA
    PatternMatch --> StabilityAnalysis : Cruciform/Triplex
    PatternMatch --> RepeatAnalysis : Slipped DNA
    PatternMatch --> iMotifScoring : i-Motif
    PatternMatch --> RLoopScoring : R-Loop
    PatternMatch --> APhilicScoring : A-philic DNA
    
    G4Hunter --> ScoreNormalization : 0-1 range
    ZSeeker --> ScoreNormalization : 0-1 range
    CurvatureAnalysis --> ScoreNormalization : 0-1 range
    StabilityAnalysis --> ScoreNormalization : 0-1 range
    RepeatAnalysis --> ScoreNormalization : 0-1 range
    iMotifScoring --> ScoreNormalization : 0-1 range
    RLoopScoring --> ScoreNormalization : 0-1 range
    APhilicScoring --> ScoreNormalization : 0-1 range
    
    ScoreNormalization --> QualityFiltering : Threshold check
    
    QualityFiltering --> PassedFilter : Score >= threshold
    QualityFiltering --> FailedFilter : Score < threshold
    
    PassedFilter --> MotifValidation : Final validation
    FailedFilter --> [*] : Discard
    
    MotifValidation --> ValidMotif : Valid
    MotifValidation --> InvalidMotif : Invalid
    
    ValidMotif --> [*] : Add to results
    InvalidMotif --> [*] : Discard
```

## ⚙️ Configuration & Class Management

```mermaid
erDiagram
    CLASSIFICATION_CONFIG {
        string class_name
        int class_id
        list subclasses
        dict scoring_methods
        dict length_limits
        dict quality_thresholds
    }
    
    CLASS_DEFINITIONS {
        string class_name
        string description
        list detection_methods
        string color_code
        bool enabled
    }
    
    MOTIF_LENGTH_LIMITS {
        int min_length
        int max_length
        string class_name
        string subclass_name
    }
    
    SCORING_METHODS {
        string algorithm_name
        string class_target
        dict parameters
        float threshold
        string normalization
    }
    
    DETECTOR_REGISTRY {
        string detector_name
        string class_target
        string module_path
        string function_name
        dict default_params
    }
    
    HYPERSCAN_PATTERNS {
        string pattern_id
        string regex_pattern
        string class_name
        int flags
        bool enabled
    }
    
    CLASSIFICATION_CONFIG ||--o{ CLASS_DEFINITIONS : "defines"
    CLASSIFICATION_CONFIG ||--o{ MOTIF_LENGTH_LIMITS : "constrains"
    CLASSIFICATION_CONFIG ||--o{ SCORING_METHODS : "specifies"
    CLASS_DEFINITIONS ||--o{ DETECTOR_REGISTRY : "maps to"
    DETECTOR_REGISTRY ||--o{ HYPERSCAN_PATTERNS : "uses"
    SCORING_METHODS ||--o{ HYPERSCAN_PATTERNS : "validates"
```

## 🚀 Performance Optimization Flow

```mermaid
flowchart TD
    START([Start Analysis]) --> CACHE_CHECK{Cache Available?}
    
    CACHE_CHECK -->|Hit| RETURN_CACHED[⚡ Return Cached Results]
    CACHE_CHECK -->|Miss| PARALLEL_SETUP[🔄 Setup Parallel Processing]
    
    PARALLEL_SETUP --> WORKER_POOL[👥 ProcessPoolExecutor<br/>Max workers = CPU cores]
    
    WORKER_POOL --> DISTRIBUTE[📦 Distribute Tasks<br/>Per-class detection]
    
    DISTRIBUTE --> HYPERSCAN_INIT[🚀 Initialize Hyperscan<br/>Pattern database compilation]
    
    HYPERSCAN_INIT --> BATCH_PROCESS[⚡ Batch Processing<br/>Vectorized operations]
    
    BATCH_PROCESS --> SIMD_SCORING[🏃‍♂️ SIMD Scoring<br/>Numba JIT compilation]
    
    SIMD_SCORING --> COLLECT_RESULTS[📊 Collect Results<br/>Merge worker outputs]
    
    COLLECT_RESULTS --> FILTER_OPTIMIZE[🎛️ Optimized Filtering<br/>Minimize overlaps]
    
    FILTER_OPTIMIZE --> CACHE_STORE[💾 Store in Cache<br/>For future requests]
    
    CACHE_STORE --> RETURN_RESULTS[✅ Return Results]
    
    RETURN_CACHED --> END([End])
    RETURN_RESULTS --> END
```

This comprehensive documentation covers all major components, workflows, and interactions within the NonBScanner system, providing clear visual representations of the tool's architecture and functionality.