# 🧬 NonBScanner Complete Flow Diagrams & Working Documentation

## 📋 Table of Contents
1. [Executive Summary](#executive-summary)
2. [System Overview](#system-overview)
3. [Architecture Mindmap](#architecture-mindmap)
4. [Complete Workflow](#complete-workflow)
5. [Technology Stack](#technology-stack)
6. [Key Features Summary](#key-features-summary)

## 📊 Executive Summary

NonBScanner (NBDFinder) is a comprehensive bioinformatics tool for detecting and analyzing Non-B DNA motifs in genomic sequences. The system combines high-performance pattern matching with scientific scoring algorithms to provide detailed analysis of structural DNA elements across 11 major classes and 22+ specialized subclasses.

### 🎯 Core Capabilities
- **Multi-Interface Access**: Web UI (Streamlit), REST API (FastAPI), and CLI
- **High-Performance Processing**: Intel Hyperscan pattern matching (40x+ speed improvement)
- **Comprehensive Detection**: 11 motif classes with 22+ subclasses
- **Advanced Visualization**: 21+ chart types and genome browser integration
- **Multiple Export Formats**: BED, BigWig, CSV, JSON, Excel

## 🏗️ System Overview

```mermaid
graph TB
    subgraph "NonBScanner System Architecture"
        subgraph "User Access Layer"
            WEB["🌐 Streamlit Web App<br/>📱 Interactive Interface<br/>📊 Real-time Visualization<br/>📤 Multi-format Export"]
            API["🔌 FastAPI REST Service<br/>🚀 High-performance API<br/>📝 Auto-documentation<br/>🔄 CORS Support"]
            CLI["💻 Command Line Interface<br/>⚡ Batch Processing<br/>📈 Progress Tracking<br/>🔧 Scripting Support"]
        end
        
        subgraph "Core Processing Engine"
            ORCHESTRATOR["🎯 Main Orchestrator<br/>🔄 Parallel Processing<br/>💾 Cache Management<br/>⚖️ Load Balancing"]
            HYPERSCAN["🚀 Hyperscan Engine<br/>⚡ Pattern Matching<br/>🧠 Intel SSE/AVX<br/>📈 40x Performance"]
            SCORING["📊 Scoring Algorithms<br/>🧬 G4Hunter<br/>🌀 Z-Seeker<br/>📐 Curvature Analysis"]
        end
        
        subgraph "Motif Detection Classes"
            CLASS1["🌊 Class 1: Curved DNA<br/>📐 A-tract Curvature<br/>🔄 Phase Analysis"]
            CLASS2["🔗 Class 2: Slipped DNA<br/>🔄 Tandem Repeats<br/>📊 STR Analysis"]
            CLASS3["✕ Class 3: Cruciform<br/>🔄 Inverted Repeats<br/>🏗️ Four-way Junction"]
            CLASS6["⭕ Class 6: G-Quadruplex<br/>🧬 7 Subclasses<br/>📊 G4Hunter Scoring"]
            CLASS7["🔵 Class 7: i-Motif<br/>🧬 C-rich Structures<br/>📊 3 Subclasses"]
            OTHERS["🧬 Classes 4,5,8,9,10,11<br/>🔄 R-Loop, Triplex<br/>🌀 Z-DNA, A-philic<br/>🔀 Hybrid, Clusters"]
        end
        
        subgraph "Output & Visualization"
            VISUALIZATION["📈 Visualization Suite<br/>📊 21+ Chart Types<br/>🎨 Interactive Plots<br/>🌐 Browser Integration"]
            EXPORT["💾 Export Engine<br/>📋 BED Format<br/>📊 BigWig Tracks<br/>📄 JSON/CSV/Excel"]
        end
    end
    
    WEB --> ORCHESTRATOR
    API --> ORCHESTRATOR
    CLI --> ORCHESTRATOR
    
    ORCHESTRATOR --> HYPERSCAN
    ORCHESTRATOR --> SCORING
    
    HYPERSCAN --> CLASS1
    HYPERSCAN --> CLASS2
    HYPERSCAN --> CLASS3
    HYPERSCAN --> CLASS6
    HYPERSCAN --> CLASS7
    HYPERSCAN --> OTHERS
    
    CLASS1 --> VISUALIZATION
    CLASS2 --> VISUALIZATION
    CLASS3 --> VISUALIZATION
    CLASS6 --> VISUALIZATION
    CLASS7 --> VISUALIZATION
    OTHERS --> VISUALIZATION
    
    VISUALIZATION --> EXPORT
```

## 🧠 Architecture Mindmap

```mermaid
mindmap
  root((NonBScanner<br/>Architecture))
    User Interfaces
      Streamlit Web App
        Interactive Dashboard
        Real-time Analysis
        Drag & Drop Upload
        Multi-format Export
      FastAPI REST Service
        RESTful Endpoints
        Auto Documentation
        JSON Responses
        CORS Support
      Command Line Interface
        Batch Processing
        Script Integration
        Progress Tracking
        Pipeline Support
    
    Core Engine
      Main Orchestrator
        Parallel Processing
        Cache Management
        Error Handling
        Resource Management
      Hyperscan Integration
        Pattern Compilation
        Fast Scanning
        Intel Optimization
        Memory Efficient
      Scoring Systems
        G4Hunter Algorithm
        Z-DNA Detection
        Curvature Analysis
        Statistical Models
    
    Motif Classes
      Structural Classes
        Curved DNA (A-tracts)
        Cruciform (Inverted repeats)
        Z-DNA (Left-handed helix)
        Triplex (Three-strand)
      Repeat Classes
        Slipped DNA (Tandems)
        R-Loop (RNA-DNA hybrid)
        Clusters (High density)
      G-rich Classes
        G-Quadruplex (7 types)
        i-Motif (C-rich)
      Complex Classes
        A-philic DNA
        Hybrid (Multi-class)
    
    Technology Stack
      Python Core
        NumPy (Arrays)
        Pandas (DataFrames)
        BioPython (Sequences)
        SciPy (Statistics)
      Performance
        Hyperscan (Patterns)
        Numba (JIT)
        Multiprocessing
        SIMD Operations
      Visualization
        Matplotlib (Static)
        Plotly (Interactive)
        Seaborn (Statistical)
        NetworkX (Graphs)
      Web Technologies
        Streamlit (UI)
        FastAPI (API)
        Uvicorn (Server)
        Pydantic (Validation)
    
    Data Flow
      Input Processing
        FASTA Parsing
        Sequence Validation
        Multi-sequence Support
        NCBI Integration
      Analysis Pipeline
        Pattern Matching
        Parallel Detection
        Score Calculation
        Quality Filtering
      Output Generation
        Result Formatting
        Statistical Analysis
        Export Conversion
        Browser Integration
```

## 🔄 Complete Workflow

```mermaid
journey
    title NonBScanner Complete Analysis Journey
    section Input & Setup
      Choose Interface: 5: User
      Upload Sequence: 4: User, System
      Configure Analysis: 4: User
      Validate Input: 5: System
    section Pattern Matching
      Initialize Hyperscan: 5: System
      Compile Patterns: 5: System
      Fast Scanning: 5: System
      Identify Candidates: 5: System
    section Parallel Detection
      Curved DNA Detection: 5: System
      Slipped DNA Detection: 5: System
      Cruciform Detection: 5: System
      R-Loop Detection: 5: System
      Triplex Detection: 5: System
      G-Quadruplex Detection: 5: System
      i-Motif Detection: 5: System
      Z-DNA Detection: 5: System
      A-philic DNA Detection: 5: System
    section Scoring & Validation
      Calculate Scores: 5: System
      Apply Algorithms: 5: System
      Normalize Results: 4: System
      Quality Filtering: 4: System
    section Post-Processing
      Hybrid Detection: 4: System
      Cluster Analysis: 4: System
      Overlap Removal: 4: System
      Statistical Analysis: 5: System
    section Visualization & Export
      Generate Plots: 5: System
      Create Tracks: 4: System
      Format Export: 4: System
      Download Results: 5: User
```

## 🛠️ Technology Stack

```mermaid
graph LR
    subgraph "Frontend Technologies"
        STREAMLIT[🌐 Streamlit<br/>Web Interface]
        HTML[🌐 HTML/CSS<br/>Styling]
        JS[⚡ JavaScript<br/>Interactivity]
    end
    
    subgraph "Backend Technologies"
        PYTHON[🐍 Python 3.8+<br/>Core Language]
        FASTAPI[⚡ FastAPI<br/>REST API Framework]
        UVICORN[🦄 Uvicorn<br/>ASGI Server]
    end
    
    subgraph "Data Processing"
        NUMPY[🔢 NumPy<br/>Numerical Computing]
        PANDAS[🐼 Pandas<br/>Data Analysis]
        SCIPY[🧮 SciPy<br/>Scientific Computing]
        BIOPYTHON[🧬 BioPython<br/>Bioinformatics]
    end
    
    subgraph "Performance Libraries"
        HYPERSCAN[🚀 Hyperscan<br/>Pattern Matching]
        NUMBA[⚡ Numba<br/>JIT Compilation]
        MULTIPROCESSING[🔄 Multiprocessing<br/>Parallel Computing]
        CONCURRENT[🔀 Concurrent.futures<br/>Thread/Process Pools]
    end
    
    subgraph "Visualization Libraries"
        MATPLOTLIB[📊 Matplotlib<br/>Static Plotting]
        PLOTLY[📈 Plotly<br/>Interactive Charts]
        SEABORN[🎨 Seaborn<br/>Statistical Plots]
        NETWORKX[🕸️ NetworkX<br/>Network Analysis]
    end
    
    subgraph "File I/O Libraries"
        OPENPYXL[📊 OpenPyXL<br/>Excel Files]
        XLSXWRITER[✍️ XlsxWriter<br/>Excel Writing]
        PATHLIB[📁 Pathlib<br/>File Handling]
    end
    
    PYTHON --> STREAMLIT
    PYTHON --> FASTAPI
    PYTHON --> NUMPY
    PYTHON --> PANDAS
    PYTHON --> HYPERSCAN
    
    FASTAPI --> UVICORN
    NUMPY --> SCIPY
    NUMPY --> NUMBA
    PANDAS --> MATPLOTLIB
    PANDAS --> PLOTLY
```

## 🌟 Key Features Summary

```mermaid
graph TB
    subgraph "Detection Capabilities"
        CLASSES["🧬 11 Motif Classes<br/>📊 22+ Subclasses<br/>🔍 Literature-validated<br/>⚡ High-performance"]
        ALGORITHMS["📊 Scoring Algorithms<br/>🧬 G4Hunter<br/>🌀 Z-Seeker<br/>📐 Curvature Analysis"]
        QUALITY["✨ Quality Control<br/>🎛️ Threshold Filtering<br/>🚫 Overlap Removal<br/>📈 Statistical Validation"]
    end
    
    subgraph "Performance Features"
        SPEED["🚀 40x Speed Improvement<br/>⚡ Intel Hyperscan<br/>🔄 Parallel Processing<br/>💾 Smart Caching"]
        SCALE["📈 Scalable Processing<br/>🔄 Batch Analysis<br/>📁 Multi-FASTA Support<br/>☁️ Memory Efficient"]
        OPTIMIZE["⚡ Performance Optimization<br/>🧠 SIMD Operations<br/>🔧 JIT Compilation<br/>📊 Resource Management"]
    end
    
    subgraph "User Experience"
        INTERFACES["🖥️ Multiple Interfaces<br/>🌐 Web Dashboard<br/>🔌 REST API<br/>💻 Command Line"]
        VISUALIZATION["📊 Rich Visualization<br/>📈 21+ Chart Types<br/>🎨 Interactive Plots<br/>🌍 Browser Integration"]
        EXPORT["💾 Multiple Export Formats<br/>📋 BED Files<br/>📊 BigWig Tracks<br/>📄 JSON/CSV/Excel"]
    end
    
    subgraph "Scientific Accuracy"
        VALIDATION["✅ Literature-based<br/>📚 Peer-reviewed Methods<br/>🔬 Biological Relevance<br/>📊 Quality Metrics"]
        STANDARDS["📏 Standard Compliance<br/>🧬 IUPAC Nomenclature<br/>📋 BED/BigWig Format<br/>🌐 Genome Browser Support"]
        REPRODUCIBLE["🔄 Reproducible Results<br/>💾 Deterministic Analysis<br/>📝 Detailed Logging<br/>🔧 Version Control"]
    end
    
    CLASSES --> SPEED
    ALGORITHMS --> SCALE
    QUALITY --> OPTIMIZE
    
    SPEED --> INTERFACES
    SCALE --> VISUALIZATION
    OPTIMIZE --> EXPORT
    
    INTERFACES --> VALIDATION
    VISUALIZATION --> STANDARDS
    EXPORT --> REPRODUCIBLE
```

## 📈 Performance Metrics

| Metric | Value | Description |
|--------|-------|-------------|
| **Speed Improvement** | 40x+ | Hyperscan vs traditional regex |
| **Motif Classes** | 11 | Major structural categories |
| **Subclasses** | 22+ | Specialized detection types |
| **Chart Types** | 21+ | Visualization options |
| **Export Formats** | 6+ | Output file types |
| **Parallel Workers** | CPU cores | Automatic scaling |
| **Memory Efficiency** | Streaming | Large file support |
| **Cache Hit Rate** | >90% | Repeat analysis speedup |

## 🔗 Documentation Links

1. **[System Architecture](SYSTEM_ARCHITECTURE.md)** - Complete system overview with main flow diagrams
2. **[Detailed Components](DETAILED_COMPONENTS.md)** - In-depth component analysis and interactions
3. **[User Workflows](USER_WORKFLOWS.md)** - User interface flows and interaction patterns
4. **[Technical Implementation](TECHNICAL_IMPLEMENTATION.md)** - Code structure and implementation details

## 📞 Quick Start Commands

```bash
# Clone repository
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner

# Install dependencies
pip install -r requirements.txt

# Launch web interface
streamlit run app.py

# Start REST API
python api.py

# Command line analysis
python cli/main.py --input sequence.fasta --output results.csv
```

## 🎯 Use Cases

- **🔬 Research**: Academic and industrial research in DNA structure
- **🧬 Genomics**: Large-scale genomic analysis pipelines
- **💊 Drug Discovery**: Targeting non-B DNA structures
- **📊 Bioinformatics**: Integration with existing analysis workflows
- **📚 Education**: Teaching structural biology and bioinformatics

This comprehensive documentation provides a complete overview of the NonBScanner tool's architecture, functionality, and usage patterns through detailed flow diagrams and technical specifications.