# 🧬 NonBScanner System Architecture & Flow Diagrams

## Overview
NonBScanner (NBDFinder) is a comprehensive Non-B DNA motif detection suite with multi-interface capabilities, high-performance pattern matching, and extensive visualization features.

## 🏗️ Main System Architecture

```mermaid
graph TB
    subgraph "User Interfaces"
        WEB[🌐 Streamlit Web App<br/>Port 8501]
        API[🔌 FastAPI REST Service<br/>Port 8000]
        CLI[💻 Command Line Interface]
    end
    
    subgraph "Core Engine"
        ORCH[🎯 Main Orchestrator<br/>all_motifs_refactored.py]
        CORE[⚡ Core Module<br/>core/]
        DETECT[🔍 Detector Registry<br/>detectors/]
    end
    
    subgraph "Processing Pipeline"
        INPUT[📄 Input Processing<br/>FASTA/Raw Text]
        HYPER[🚀 Hyperscan Engine<br/>Pattern Matching]
        SCORE[📊 Scoring Algorithms<br/>G4Hunter, Z-Seeker, etc.]
        FILTER[🎛️ Post-Processing<br/>Overlap removal, clustering]
    end
    
    subgraph "Motif Classes (10 Main + 22 Subclasses)"
        C1[Class 1: Curved DNA]
        C2[Class 2: Slipped DNA]
        C3[Class 3: Cruciform]
        C4[Class 4: R-Loop]
        C5[Class 5: Triplex]
        C6[Class 6: G-Quadruplex]
        C7[Class 7: i-Motif]
        C8[Class 8: Z-DNA]
        C9[Class 9: A-philic DNA]
        C10[Class 10: Hybrid]
        C11[Class 11: Clusters]
    end
    
    subgraph "Output & Visualization"
        VIZ[📈 Visualization Suite<br/>viz/]
        EXPORT[💾 Export Formats<br/>BED, BigWig, CSV, JSON]
        BROWSER[🌍 Genome Browser<br/>UCSC, IGV, JBrowse]
    end
    
    WEB --> ORCH
    API --> ORCH
    CLI --> ORCH
    
    ORCH --> CORE
    ORCH --> DETECT
    
    CORE --> INPUT
    INPUT --> HYPER
    HYPER --> SCORE
    SCORE --> FILTER
    
    FILTER --> C1
    FILTER --> C2
    FILTER --> C3
    FILTER --> C4
    FILTER --> C5
    FILTER --> C6
    FILTER --> C7
    FILTER --> C8
    FILTER --> C9
    FILTER --> C10
    FILTER --> C11
    
    C1 --> VIZ
    C2 --> VIZ
    C3 --> VIZ
    C4 --> VIZ
    C5 --> VIZ
    C6 --> VIZ
    C7 --> VIZ
    C8 --> VIZ
    C9 --> VIZ
    C10 --> VIZ
    C11 --> VIZ
    
    VIZ --> EXPORT
    VIZ --> BROWSER
```

## 🔄 Core Processing Workflow

```mermaid
flowchart TD
    START([🚀 Start Analysis]) --> INPUT{📥 Input Type?}
    
    INPUT -->|FASTA File| PARSE1[📄 Parse FASTA<br/>utils.parse_fasta]
    INPUT -->|Raw Text| PARSE2[✏️ Process Raw Text<br/>Sequence validation]
    INPUT -->|NCBI ID| FETCH[🌐 NCBI Fetch<br/>Bio.Entrez]
    
    PARSE1 --> VALIDATE[✅ Sequence Validation<br/>ATGC check, length > 10bp]
    PARSE2 --> VALIDATE
    FETCH --> VALIDATE
    
    VALIDATE -->|Invalid| ERROR[❌ Error: Invalid Sequence]
    VALIDATE -->|Valid| CACHE{💾 Check Cache?}
    
    CACHE -->|Hit| CACHED[⚡ Return Cached Results]
    CACHE -->|Miss| PARALLEL[🔄 Parallel Detection<br/>ProcessPoolExecutor]
    
    PARALLEL --> MOTIF1[🧬 Class 1: Curved DNA<br/>A-tract detection]
    PARALLEL --> MOTIF2[🔗 Class 2: Slipped DNA<br/>Tandem repeats]
    PARALLEL --> MOTIF3[✕ Class 3: Cruciform<br/>Inverted repeats]
    PARALLEL --> MOTIF4[🔄 Class 4: R-Loop<br/>RNA-DNA hybrids]
    PARALLEL --> MOTIF5[🔺 Class 5: Triplex<br/>Three-strand formation]
    PARALLEL --> MOTIF6[⭕ Class 6: G-Quadruplex<br/>G4Hunter scoring]
    PARALLEL --> MOTIF7[🔵 Class 7: i-Motif<br/>C-rich structures]
    PARALLEL --> MOTIF8[🔄 Class 8: Z-DNA<br/>Left-handed helix]
    PARALLEL --> MOTIF9[🧬 Class 9: A-philic DNA<br/>A-tract affinity]
    
    MOTIF1 --> COLLECT[📊 Collect All Results]
    MOTIF2 --> COLLECT
    MOTIF3 --> COLLECT
    MOTIF4 --> COLLECT
    MOTIF5 --> COLLECT
    MOTIF6 --> COLLECT
    MOTIF7 --> COLLECT
    MOTIF8 --> COLLECT
    MOTIF9 --> COLLECT
    
    COLLECT --> HYBRID[🔀 Class 10: Hybrid Detection<br/>Multi-class overlaps]
    HYBRID --> CLUSTER[🌟 Class 11: Cluster Detection<br/>High-density regions]
    
    CLUSTER --> FILTER{🎛️ Apply Filters?}
    FILTER -->|nonoverlap=True| OVERLAP[🚫 Remove Overlaps<br/>Per-class best selection]
    FILTER -->|nonoverlap=False| STANDARD[📋 Standardize Output]
    
    OVERLAP --> STANDARD
    STANDARD --> STATS[📈 Calculate Statistics<br/>GC%, scores, lengths]
    STATS --> RETURN[✅ Return Results]
    
    CACHED --> RETURN
    ERROR --> END([🏁 End])
    RETURN --> END
```

## 🏛️ Module Dependencies

```mermaid
graph LR
    subgraph "Frontend Layer"
        APP[app.py<br/>Streamlit UI]
        API_MOD[api.py<br/>FastAPI]
        CLI_MOD[cli/main.py<br/>CLI Interface]
    end
    
    subgraph "Orchestration Layer"
        MAIN_ORCH[all_motifs_refactored.py<br/>Main Orchestrator]
        ORCH_MOD[orchestrators/<br/>Specialized Orchestrators]
    end
    
    subgraph "Core Processing"
        CORE_MOD[core/<br/>Core Engine]
        DETECT_MOD[detectors/<br/>Motif Detectors]
        MOTIFS_MOD[motifs/<br/>Detection Algorithms]
    end
    
    subgraph "Utilities & I/O"
        UTILS[utils.py<br/>Basic Utilities]
        NBDIO[nbdio/<br/>I/O Handling]
        CONFIG[classification_config.py<br/>Motif Configuration]
    end
    
    subgraph "Visualization"
        VIZ_MOD[viz/<br/>Plotting & Browser]
    end
    
    subgraph "External Libraries"
        STREAMLIT[streamlit]
        FASTAPI[fastapi]
        HYPERSCAN[hyperscan]
        BIOPYTHON[biopython]
        NUMPY[numpy/pandas]
        PLOTLY[plotly/matplotlib]
    end
    
    APP --> MAIN_ORCH
    API_MOD --> MAIN_ORCH
    CLI_MOD --> MAIN_ORCH
    
    MAIN_ORCH --> CORE_MOD
    MAIN_ORCH --> DETECT_MOD
    ORCH_MOD --> CORE_MOD
    
    CORE_MOD --> MOTIFS_MOD
    DETECT_MOD --> MOTIFS_MOD
    
    CORE_MOD --> UTILS
    MOTIFS_MOD --> UTILS
    MAIN_ORCH --> CONFIG
    
    APP --> VIZ_MOD
    VIZ_MOD --> PLOTLY
    
    MOTIFS_MOD --> NBDIO
    NBDIO --> BIOPYTHON
    
    APP --> STREAMLIT
    API_MOD --> FASTAPI
    CORE_MOD --> HYPERSCAN
    CORE_MOD --> NUMPY
```

## 🎯 Motif Detection Pipeline

```mermaid
sequenceDiagram
    participant User
    participant Interface as Web/API/CLI
    participant Orchestrator
    participant Hyperscan
    participant Detector as Motif Detector
    participant Scorer as Scoring Engine
    participant PostProc as Post-Processor
    
    User->>Interface: Submit DNA Sequence
    Interface->>Orchestrator: Process Request
    
    Orchestrator->>Orchestrator: Validate Sequence
    Orchestrator->>Orchestrator: Check Cache
    
    alt Cache Miss
        Orchestrator->>Hyperscan: Initialize Pattern Database
        Hyperscan->>Detector: Fast Pattern Matching
        
        par Parallel Detection
            Detector->>Scorer: Class 1: Curved DNA
            Detector->>Scorer: Class 2: Slipped DNA
            Detector->>Scorer: Class 3: Cruciform
            Detector->>Scorer: Class 4: R-Loop
            Detector->>Scorer: Class 5: Triplex
            Detector->>Scorer: Class 6: G-Quadruplex
            Detector->>Scorer: Class 7: i-Motif
            Detector->>Scorer: Class 8: Z-DNA
            Detector->>Scorer: Class 9: A-philic DNA
        end
        
        Scorer->>PostProc: Raw Motif Results
        PostProc->>PostProc: Hybrid Detection (Class 10)
        PostProc->>PostProc: Cluster Detection (Class 11)
        PostProc->>PostProc: Overlap Removal
        PostProc->>PostProc: Quality Filtering
        PostProc->>Orchestrator: Processed Results
        Orchestrator->>Orchestrator: Cache Results
    else Cache Hit
        Orchestrator->>Orchestrator: Retrieve Cached Results
    end
    
    Orchestrator->>Interface: Return Formatted Results
    Interface->>User: Display/Export Results
```

## 📊 Data Flow Architecture

```mermaid
graph TD
    subgraph "Input Sources"
        FILE[📁 FASTA Files]
        TEXT[✏️ Raw Text Input]
        NCBI[🌐 NCBI Database]
        EXAMPLE[🧪 Example Data]
    end
    
    subgraph "Input Processing"
        PARSER[📝 FASTA Parser<br/>utils.parse_fasta]
        VALIDATOR[✅ Sequence Validator<br/>ATGC check, min length]
        NORMALIZER[🔄 Sequence Normalizer<br/>Uppercase, clean]
    end
    
    subgraph "Pattern Matching Engine"
        REGEX[🔍 Regex Registry<br/>core/regex_registry]
        HYPERSCAN_DB[🚀 Hyperscan Database<br/>Compiled patterns]
        MATCHER[⚡ Pattern Matcher<br/>High-speed scanning]
    end
    
    subgraph "Scoring Systems"
        G4HUNTER[🧬 G4Hunter<br/>G-quadruplex scoring]
        ZSEEKER[🔄 Z-Seeker<br/>Z-DNA detection]
        CURVATURE[📐 Curvature Analysis<br/>A-tract bending]
        STABILITY[🔬 Stability Scoring<br/>Thermodynamic models]
    end
    
    subgraph "Classification & Filtering"
        CLASSIFIER[🏷️ Motif Classifier<br/>22+ subclasses]
        OVERLAP[🚫 Overlap Filter<br/>Per-class optimization]
        QUALITY[✨ Quality Filter<br/>Score/length thresholds]
        MERGER[🔗 Nearby Merger<br/>Close motif combination]
    end
    
    subgraph "Output Generation"
        FORMATTER[📋 Result Formatter<br/>Standardized output]
        STATISTICS[📊 Statistics Calculator<br/>Summary metrics]
        EXPORTER[💾 Multi-format Export<br/>BED, CSV, JSON, BigWig]
        VISUALIZER[📈 Visualization<br/>Charts, plots, tracks]
    end
    
    FILE --> PARSER
    TEXT --> PARSER
    NCBI --> PARSER
    EXAMPLE --> PARSER
    
    PARSER --> VALIDATOR
    VALIDATOR --> NORMALIZER
    
    NORMALIZER --> REGEX
    REGEX --> HYPERSCAN_DB
    HYPERSCAN_DB --> MATCHER
    
    MATCHER --> G4HUNTER
    MATCHER --> ZSEEKER
    MATCHER --> CURVATURE
    MATCHER --> STABILITY
    
    G4HUNTER --> CLASSIFIER
    ZSEEKER --> CLASSIFIER
    CURVATURE --> CLASSIFIER
    STABILITY --> CLASSIFIER
    
    CLASSIFIER --> OVERLAP
    OVERLAP --> QUALITY
    QUALITY --> MERGER
    
    MERGER --> FORMATTER
    FORMATTER --> STATISTICS
    STATISTICS --> EXPORTER
    STATISTICS --> VISUALIZER
```

## 🌐 API Interaction Flow

```mermaid
graph TB
    subgraph "Client Applications"
        WEBCLIENT[🌐 Web Browser<br/>Streamlit Interface]
        APICLIENT[🔌 API Client<br/>REST requests]
        CURL[💻 cURL/Postman<br/>Direct API calls]
        PYTHON[🐍 Python Scripts<br/>requests library]
    end
    
    subgraph "API Endpoints"
        HEALTH[📊 /api/v1/health<br/>System status]
        CLASSES[📋 /api/v1/classes<br/>Motif class info]
        ANALYZE[🔬 /api/v1/analyze<br/>Full analysis]
        CLASS_ANALYZE[🎯 /api/v1/analyze/{class}<br/>Specific class]
        STATS[📈 /api/v1/stats<br/>Usage statistics]
        DOCS[📚 /docs<br/>API documentation]
    end
    
    subgraph "Middleware & Processing"
        CORS[🔗 CORS Middleware<br/>Cross-origin support]
        AUTH[🔐 Authentication<br/>Optional security]
        RATE[⏱️ Rate Limiting<br/>Traffic control]
        CACHE[💾 Response Cache<br/>Performance boost]
    end
    
    subgraph "Core Processing"
        VALIDATOR_API[✅ Request Validator<br/>Pydantic models]
        ORCHESTRATOR_API[🎯 Analysis Engine<br/>all_motifs_refactored]
        FORMATTER_API[📋 Response Formatter<br/>JSON standardization]
    end
    
    WEBCLIENT --> ANALYZE
    APICLIENT --> ANALYZE
    CURL --> HEALTH
    PYTHON --> CLASS_ANALYZE
    
    HEALTH --> CORS
    CLASSES --> CORS
    ANALYZE --> CORS
    CLASS_ANALYZE --> CORS
    STATS --> CORS
    DOCS --> CORS
    
    CORS --> AUTH
    AUTH --> RATE
    RATE --> CACHE
    
    CACHE --> VALIDATOR_API
    VALIDATOR_API --> ORCHESTRATOR_API
    ORCHESTRATOR_API --> FORMATTER_API
    
    FORMATTER_API --> CACHE
```