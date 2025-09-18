# 🖥️ NonBScanner User Interface & Workflow Diagrams

## 📱 User Interface Overview

```mermaid
graph TB
    subgraph "User Entry Points"
        WEB_BROWSER[🌐 Web Browser<br/>http://localhost:8501]
        API_CLIENT[🔌 API Client<br/>http://localhost:8000]
        CLI_TERMINAL[💻 Command Line<br/>python cli/main.py]
    end
    
    subgraph "Streamlit Web Interface"
        HOME_PAGE[🏠 Home Page<br/>Overview & Documentation]
        UPLOAD_PAGE[📤 Upload & Analyze<br/>Sequence Input]
        RESULTS_PAGE[📊 Results<br/>Analysis Output]
        DOWNLOAD_PAGE[💾 Download<br/>Export Options]
        DOCS_PAGE[📚 Documentation<br/>Scientific References]
    end
    
    subgraph "Input Methods"
        FILE_UPLOAD[📁 FASTA File Upload<br/>Drag & drop interface]
        TEXT_INPUT[✏️ Paste Sequence<br/>Raw text input]
        EXAMPLE_DATA[🧪 Example Data<br/>Pre-loaded samples]
        NCBI_FETCH[🌐 NCBI Fetch<br/>Database retrieval]
    end
    
    subgraph "Analysis Options"
        CLASS_SELECT[🎯 Class Selection<br/>Choose motif classes]
        PARAM_CONFIG[⚙️ Parameter Configuration<br/>Scoring thresholds]
        OVERLAP_FILTER[🎛️ Overlap Filtering<br/>nonoverlap option]
        HOTSPOT_DETECT[🌟 Hotspot Detection<br/>Cluster analysis]
    end
    
    subgraph "Results Visualization"
        SUMMARY_STATS[📈 Summary Statistics<br/>Counts, distributions]
        MOTIF_TABLE[📋 Motif Table<br/>Detailed results]
        INTERACTIVE_PLOTS[📊 Interactive Plots<br/>21+ chart types]
        GENOME_TRACKS[🧬 Genome Tracks<br/>Browser-style view]
    end
    
    subgraph "Export Options"
        CSV_EXPORT[📊 CSV Export<br/>Tabular data]
        BED_EXPORT[📋 BED Format<br/>Genome annotations]
        BIGWIG_EXPORT[📈 BigWig Format<br/>Signal tracks]
        JSON_EXPORT[📄 JSON Export<br/>Structured data]
        EXCEL_EXPORT[📊 Excel Export<br/>Multi-sheet workbook]
        PNG_EXPORT[🖼️ PNG/PDF Export<br/>Publication plots]
    end
    
    WEB_BROWSER --> HOME_PAGE
    HOME_PAGE --> UPLOAD_PAGE
    HOME_PAGE --> DOCS_PAGE
    
    UPLOAD_PAGE --> FILE_UPLOAD
    UPLOAD_PAGE --> TEXT_INPUT
    UPLOAD_PAGE --> EXAMPLE_DATA
    UPLOAD_PAGE --> NCBI_FETCH
    
    FILE_UPLOAD --> CLASS_SELECT
    TEXT_INPUT --> CLASS_SELECT
    EXAMPLE_DATA --> CLASS_SELECT
    NCBI_FETCH --> CLASS_SELECT
    
    CLASS_SELECT --> PARAM_CONFIG
    PARAM_CONFIG --> OVERLAP_FILTER
    OVERLAP_FILTER --> HOTSPOT_DETECT
    
    HOTSPOT_DETECT --> RESULTS_PAGE
    
    RESULTS_PAGE --> SUMMARY_STATS
    RESULTS_PAGE --> MOTIF_TABLE
    RESULTS_PAGE --> INTERACTIVE_PLOTS
    RESULTS_PAGE --> GENOME_TRACKS
    
    RESULTS_PAGE --> DOWNLOAD_PAGE
    
    DOWNLOAD_PAGE --> CSV_EXPORT
    DOWNLOAD_PAGE --> BED_EXPORT
    DOWNLOAD_PAGE --> BIGWIG_EXPORT
    DOWNLOAD_PAGE --> JSON_EXPORT
    DOWNLOAD_PAGE --> EXCEL_EXPORT
    DOWNLOAD_PAGE --> PNG_EXPORT
    
    API_CLIENT --> API_ENDPOINTS[🔌 REST API Endpoints]
    CLI_TERMINAL --> CLI_COMMANDS[💻 Command Line Tools]
```

## 🔄 Complete User Workflow

```mermaid
journey
    title NonBScanner User Journey
    section Sequence Input
      Open Web Interface: 5: User
      Choose Input Method: 4: User
      Upload FASTA File: 3: User, System
      Validate Sequence: 4: System
    section Analysis Configuration
      Select Motif Classes: 4: User
      Configure Parameters: 3: User
      Set Filtering Options: 3: User
      Start Analysis: 5: User, System
    section Processing
      Initialize Hyperscan: 5: System
      Parallel Detection: 5: System
      Score Calculation: 5: System
      Post-processing: 4: System
    section Results Review
      View Summary Stats: 5: User
      Explore Motif Table: 4: User
      Interactive Plotting: 5: User
      Genome Track View: 4: User
    section Export & Download
      Choose Export Format: 4: User
      Generate Files: 4: System
      Download Results: 5: User
```

## 🌐 API Usage Patterns

```mermaid
sequenceDiagram
    participant Client
    participant API as FastAPI Server
    participant Cache
    participant Engine as Analysis Engine
    
    Note over Client,Engine: Health Check & System Info
    Client->>API: GET /api/v1/health
    API-->>Client: System status
    
    Client->>API: GET /api/v1/classes
    API-->>Client: Available motif classes
    
    Note over Client,Engine: Full Sequence Analysis
    Client->>API: POST /api/v1/analyze
    Note right of Client: {<br/>  "sequence": "ATCG...",<br/>  "sequence_name": "test",<br/>  "nonoverlap": true<br/>}
    
    API->>Cache: Check cache
    alt Cache Hit
        Cache-->>API: Cached results
    else Cache Miss
        API->>Engine: Process sequence
        Engine-->>API: Analysis results
        API->>Cache: Store results
    end
    
    API-->>Client: Complete analysis
    Note left of API: {<br/>  "motifs": [...],<br/>  "statistics": {...},<br/>  "metadata": {...}<br/>}
    
    Note over Client,Engine: Class-Specific Analysis
    Client->>API: POST /api/v1/analyze/6
    Note right of Client: G-Quadruplex only
    
    API->>Engine: Class-specific detection
    Engine-->>API: G4 results only
    API-->>Client: Class 6 motifs
    
    Note over Client,Engine: Usage Statistics
    Client->>API: GET /api/v1/stats
    API-->>Client: API usage metrics
```

## 💻 Command Line Interface Flow

```mermaid
flowchart TD
    START([CLI Start]) --> ARGS_PARSE[📝 Parse Arguments<br/>argparse]
    
    ARGS_PARSE --> INPUT_CHECK{📥 Input Type?}
    
    INPUT_CHECK -->|File| FILE_READ[📁 Read FASTA File<br/>nbdio.fasta.parse]
    INPUT_CHECK -->|Directory| DIR_SCAN[📂 Scan Directory<br/>Batch processing]
    INPUT_CHECK -->|Stdin| STDIN_READ[⌨️ Read from stdin<br/>Pipe support]
    
    FILE_READ --> VALIDATE_SEQ[✅ Validate Sequences]
    DIR_SCAN --> VALIDATE_SEQ
    STDIN_READ --> VALIDATE_SEQ
    
    VALIDATE_SEQ --> CONFIG_LOAD[⚙️ Load Configuration<br/>CLI parameters]
    
    CONFIG_LOAD --> PROGRESS_INIT[📊 Initialize Progress<br/>tqdm progress bar]
    
    PROGRESS_INIT --> BATCH_PROCESS[🔄 Batch Processing<br/>Multiple sequences]
    
    BATCH_PROCESS --> ANALYZE_SEQ[🧬 Analyze Sequence<br/>all_motifs_refactored]
    
    ANALYZE_SEQ --> UPDATE_PROGRESS[📈 Update Progress<br/>Real-time feedback]
    
    UPDATE_PROGRESS --> MORE_SEQS{More Sequences?}
    
    MORE_SEQS -->|Yes| ANALYZE_SEQ
    MORE_SEQS -->|No| OUTPUT_FORMAT[📋 Format Output<br/>User-specified format]
    
    OUTPUT_FORMAT --> WRITE_RESULTS[💾 Write Results<br/>File or stdout]
    
    WRITE_RESULTS --> SUMMARY[📊 Print Summary<br/>Final statistics]
    
    SUMMARY --> END([CLI End])
```

## 📊 Visualization Interface Flow

```mermaid
stateDiagram-v2
    [*] --> AnalysisComplete
    
    AnalysisComplete --> PlotSelection : User chooses visualization
    
    PlotSelection --> DistributionPlots : Motif counts/lengths
    PlotSelection --> PositionPlots : Track-style plots
    PlotSelection --> ScorePlots : Score distributions
    PlotSelection --> InteractivePlots : Plotly dashboards
    PlotSelection --> StatisticalPlots : Advanced analytics
    
    DistributionPlots --> BarChart : Class counts
    DistributionPlots --> PieChart : Proportions
    DistributionPlots --> Histogram : Length distribution
    
    PositionPlots --> TrackPlot : Genome browser style
    PositionPlots --> DensityPlot : Motif density
    PositionPlots --> ManhattanPlot : Position vs score
    
    ScorePlots --> BoxPlot : Score by class
    ScorePlots --> ViolinPlot : Score distributions
    ScorePlots --> ScatterPlot : Score correlations
    
    InteractivePlots --> SunburstChart : Hierarchical data
    InteractivePlots --> TreemapChart : Proportional data
    InteractivePlots --> NetworkGraph : Motif interactions
    
    StatisticalPlots --> HeatmapPlot : Correlation matrix
    StatisticalPlots --> TSNEPlot : Dimensionality reduction
    StatisticalPlots --> VennDiagram : Class overlaps
    
    BarChart --> ExportOptions
    PieChart --> ExportOptions
    Histogram --> ExportOptions
    TrackPlot --> ExportOptions
    DensityPlot --> ExportOptions
    ManhattanPlot --> ExportOptions
    BoxPlot --> ExportOptions
    ViolinPlot --> ExportOptions
    ScatterPlot --> ExportOptions
    SunburstChart --> ExportOptions
    TreemapChart --> ExportOptions
    NetworkGraph --> ExportOptions
    HeatmapPlot --> ExportOptions
    TSNEPlot --> ExportOptions
    VennDiagram --> ExportOptions
    
    ExportOptions --> SavePNG : Static image
    ExportOptions --> SaveSVG : Vector format
    ExportOptions --> SaveHTML : Interactive plot
    ExportOptions --> SavePDF : Publication quality
    
    SavePNG --> [*]
    SaveSVG --> [*]
    SaveHTML --> [*]
    SavePDF --> [*]
```

## 🎛️ Configuration & Parameter Flow

```mermaid
graph LR
    subgraph "Default Configuration"
        DEFAULT_CLASSES[🎯 Default Classes<br/>All 11 classes enabled]
        DEFAULT_PARAMS[⚙️ Default Parameters<br/>Standard thresholds]
        DEFAULT_SCORING[📊 Default Scoring<br/>Literature-based]
    end
    
    subgraph "User Customization"
        CLASS_SELECTION[☑️ Class Selection<br/>Enable/disable classes]
        PARAM_OVERRIDE[🔧 Parameter Override<br/>Custom thresholds]
        SCORING_ADJUST[📈 Scoring Adjustment<br/>Weight modifications]
    end
    
    subgraph "Runtime Configuration"
        MERGE_CONFIG[🔀 Merge Configuration<br/>User + defaults]
        VALIDATE_CONFIG[✅ Validate Configuration<br/>Range checks]
        APPLY_CONFIG[⚡ Apply Configuration<br/>Analysis pipeline]
    end
    
    subgraph "Advanced Options"
        OVERLAP_HANDLING[🎛️ Overlap Handling<br/>nonoverlap setting]
        CACHE_SETTINGS[💾 Cache Settings<br/>Enable/disable caching]
        PARALLEL_CONFIG[🔄 Parallel Config<br/>Worker count]
        OUTPUT_FORMAT[📋 Output Format<br/>Export preferences]
    end
    
    DEFAULT_CLASSES --> MERGE_CONFIG
    DEFAULT_PARAMS --> MERGE_CONFIG
    DEFAULT_SCORING --> MERGE_CONFIG
    
    CLASS_SELECTION --> MERGE_CONFIG
    PARAM_OVERRIDE --> MERGE_CONFIG
    SCORING_ADJUST --> MERGE_CONFIG
    
    MERGE_CONFIG --> VALIDATE_CONFIG
    VALIDATE_CONFIG --> APPLY_CONFIG
    
    OVERLAP_HANDLING --> APPLY_CONFIG
    CACHE_SETTINGS --> APPLY_CONFIG
    PARALLEL_CONFIG --> APPLY_CONFIG
    OUTPUT_FORMAT --> APPLY_CONFIG
```

## 🚀 Performance Monitoring Flow

```mermaid
sequenceDiagram
    participant User
    participant Interface
    participant Monitor as Performance Monitor
    participant Cache
    participant Engine
    
    User->>Interface: Start Analysis
    Interface->>Monitor: Initialize tracking
    
    Monitor->>Monitor: Start timer
    Monitor->>Monitor: Track memory usage
    
    Interface->>Cache: Check cache
    Cache-->>Monitor: Cache hit/miss
    
    alt Cache Miss
        Interface->>Engine: Process sequence
        Engine-->>Monitor: Processing status
        Monitor->>Monitor: Update progress
        Engine-->>Interface: Results
    else Cache Hit
        Cache-->>Interface: Cached results
        Monitor->>Monitor: Record cache efficiency
    end
    
    Interface->>Monitor: Analysis complete
    Monitor->>Monitor: Calculate metrics
    
    Monitor->>Interface: Performance data
    Interface->>User: Results + metrics
    
    Note over User,Engine: Performance Metrics:<br/>- Processing time<br/>- Memory usage<br/>- Cache efficiency<br/>- Motif detection rate
```

This comprehensive documentation provides detailed visual representations of all user interfaces, workflows, and interaction patterns within the NonBScanner system.