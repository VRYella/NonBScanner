# 📚 NonBScanner Documentation Index

Welcome to the comprehensive documentation for NonBScanner (NBDFinder) - a state-of-the-art Non-B DNA motif detection suite.

## 📋 Documentation Structure

### 🏗️ [System Architecture](SYSTEM_ARCHITECTURE.md)
**Main system overview with core flow diagrams**
- Complete system architecture diagram
- Core processing workflow
- Module dependencies 
- Motif detection pipeline
- Data flow architecture
- API interaction flow

### 🔧 [Detailed Components](DETAILED_COMPONENTS.md) 
**In-depth component analysis and technical details**
- Motif detection classes hierarchy (mindmap)
- Hyperscan pattern matching flow
- Visualization pipeline details
- Function & library dependencies
- Scoring algorithm workflows
- Configuration & class management

### 🖥️ [User Workflows](USER_WORKFLOWS.md)
**User interface flows and interaction patterns**
- Complete user interface overview
- User journey mapping
- API usage patterns
- Command line interface flows
- Visualization interface flows
- Configuration & parameter management

### 🛠️ [Technical Implementation](TECHNICAL_IMPLEMENTATION.md)
**Code structure and implementation details**
- File structure & script organization
- Key function & class architecture
- Performance optimization flows
- Data structure flows
- Configuration system architecture
- Core algorithm implementations

### 📖 [Complete Overview](README.md)
**Executive summary and comprehensive documentation**
- Executive summary
- Complete system overview
- Architecture mindmap
- Technology stack overview
- Key features summary
- Performance metrics

## 🚀 Quick Navigation

| Topic | File | Description |
|-------|------|-------------|
| **🏗️ Architecture** | [SYSTEM_ARCHITECTURE.md](SYSTEM_ARCHITECTURE.md) | Main system flow diagrams |
| **🔍 Components** | [DETAILED_COMPONENTS.md](DETAILED_COMPONENTS.md) | Component analysis & interactions |
| **👥 User Flows** | [USER_WORKFLOWS.md](USER_WORKFLOWS.md) | Interface workflows & patterns |
| **💻 Implementation** | [TECHNICAL_IMPLEMENTATION.md](TECHNICAL_IMPLEMENTATION.md) | Code structure & algorithms |
| **📋 Overview** | [README.md](README.md) | Complete system summary |

## 🎯 What You'll Find

### 🧬 **Motif Detection Coverage**
- **11 Major Classes**: Curved DNA, Slipped DNA, Cruciform, R-Loop, Triplex, G-Quadruplex, i-Motif, Z-DNA, A-philic DNA, Hybrid, Clusters
- **22+ Subclasses**: Detailed specializations for each major class
- **Literature-validated**: Peer-reviewed detection algorithms

### ⚡ **Performance Features**  
- **40x Speed**: Intel Hyperscan pattern matching engine
- **Parallel Processing**: Multi-core CPU utilization
- **Smart Caching**: Efficient result reuse
- **Memory Optimization**: Streaming analysis for large files

### 🌐 **Multiple Interfaces**
- **Web Interface**: Streamlit dashboard (Port 8501)
- **REST API**: FastAPI service (Port 8000)  
- **Command Line**: Batch processing CLI

### 📊 **Rich Visualization**
- **21+ Chart Types**: From basic distributions to advanced t-SNE plots
- **Interactive Plots**: Plotly-powered dashboards
- **Genome Browser**: UCSC, IGV, JBrowse integration
- **Export Options**: PNG, SVG, PDF, HTML formats

### 💾 **Multiple Export Formats**
- **BED Files**: Genome coordinate annotations
- **BigWig Tracks**: Signal visualization
- **CSV/Excel**: Tabular data analysis
- **JSON**: Structured data exchange

## 🔧 Quick Start

```bash
# Launch web interface
streamlit run app.py

# Start REST API  
python api.py

# Command line analysis
python cli/main.py --input sequence.fasta --output results.csv
```

## 📈 System Capabilities

| Feature | Capability | Performance |
|---------|------------|-------------|
| **Detection Classes** | 11 major + 22+ subclasses | Comprehensive coverage |
| **Speed Improvement** | 40x+ vs traditional methods | Intel Hyperscan optimization |
| **Parallel Processing** | CPU core scaling | Automatic load balancing |
| **Visualization Types** | 21+ chart varieties | Interactive & static |
| **Export Formats** | 6+ file types | Industry standards |
| **Memory Usage** | Streaming analysis | Large file support |

## 📚 Understanding the Diagrams

The documentation uses **Mermaid diagrams** throughout for clear visual representation:

- **📊 Flowcharts**: Process flows and decision trees
- **🏗️ Architecture Diagrams**: System component relationships  
- **🧠 Mindmaps**: Hierarchical information organization
- **📈 Sequence Diagrams**: Interaction timelines
- **🗃️ Entity Diagrams**: Data structure relationships
- **🚀 State Diagrams**: Process state transitions

## 🎯 Target Audiences

### 🔬 **Researchers**
- Understand detection algorithms and scientific accuracy
- Review visualization capabilities for publication
- Explore API integration for pipelines

### 💻 **Developers** 
- Study system architecture and code organization
- Learn performance optimization techniques
- Understand module dependencies and interfaces

### 👥 **Users**
- Navigate interface workflows and features
- Understand input/output formats and options
- Learn analysis configuration and customization

### 🏫 **Educators**
- Use diagrams for teaching bioinformatics concepts
- Demonstrate real-world software architecture
- Show integration of multiple technologies

---

*This documentation provides comprehensive coverage of the NonBScanner tool through detailed flow diagrams, technical specifications, and visual representations of all system components and workflows.*