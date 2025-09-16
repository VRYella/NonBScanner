#!/bin/bash
# NBDFinder Services Launcher
# ===========================
# 
# This script launches the NBDFinder web services:
# - Streamlit web interface on port 8501
# - FastAPI REST API on port 8000

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 NBDFinder Services Launcher${NC}"
echo -e "${BLUE}=================================${NC}"
echo ""

# Check if Python packages are installed
echo -e "${YELLOW}📦 Checking dependencies...${NC}"
python3 -c "import streamlit, fastapi, uvicorn" 2>/dev/null || {
    echo -e "${RED}❌ Missing dependencies. Installing...${NC}"
    pip install -r requirements.txt
}
echo -e "${GREEN}✅ Dependencies OK${NC}"

# Function to start services
start_services() {
    echo ""
    echo -e "${YELLOW}🚀 Starting NBDFinder services...${NC}"
    
    # Start FastAPI in background
    echo -e "${BLUE}Starting REST API on http://localhost:8000${NC}"
    python3 api.py &
    API_PID=$!
    
    # Wait a moment for API to start
    sleep 3
    
    # Start Streamlit in background
    echo -e "${BLUE}Starting Web Interface on http://localhost:8501${NC}"
    streamlit run app.py --server.port 8501 --server.headless true &
    STREAMLIT_PID=$!
    
    # Wait a moment for services to start
    sleep 5
    
    echo ""
    echo -e "${GREEN}🎉 NBDFinder services are running!${NC}"
    echo ""
    echo -e "${BLUE}📱 Web Interface:${NC} http://localhost:8501"
    echo -e "${BLUE}🚀 REST API:${NC} http://localhost:8000"
    echo -e "${BLUE}📚 API Docs:${NC} http://localhost:8000/docs"
    echo ""
    echo -e "${YELLOW}📋 Available Services:${NC}"
    echo "  • Interactive web interface with visualization suite"
    echo "  • REST API for programmatic access"
    echo "  • 10 major Non-B DNA classes"
    echo "  • 22+ specialized subclasses"
    echo "  • Hyperscan-accelerated detection"
    echo "  • Comprehensive analysis and export tools"
    echo ""
    echo -e "${YELLOW}💡 Quick Examples:${NC}"
    echo "  # Test API health"
    echo "  curl http://localhost:8000/api/v1/health"
    echo ""
    echo "  # Get motif class info"
    echo "  curl http://localhost:8000/api/v1/motif-info"
    echo ""
    echo "  # Analyze sequence"
    echo '  curl -X POST http://localhost:8000/api/v1/analyze \\'
    echo '    -H "Content-Type: application/json" \\'
    echo '    -d '"'"'{"sequence": "GGGTTAGGGTTAGGGTTAGGG", "sequence_name": "test"}'"'"''
    echo ""
    echo -e "${YELLOW}🛑 Press Ctrl+C to stop all services${NC}"
    
    # Handle cleanup on exit
    cleanup() {
        echo ""
        echo -e "${YELLOW}🛑 Stopping services...${NC}"
        kill $API_PID $STREAMLIT_PID 2>/dev/null || true
        echo -e "${GREEN}✅ Services stopped${NC}"
        exit 0
    }
    
    trap cleanup SIGINT SIGTERM
    
    # Keep script running
    wait
}

# Function to show help
show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  start         Start both web interface and API services (default)"
    echo "  api-only      Start only the REST API service"
    echo "  web-only      Start only the Streamlit web interface"
    echo "  help          Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0              # Start both services"
    echo "  $0 start        # Start both services"
    echo "  $0 api-only     # Start only REST API"
    echo "  $0 web-only     # Start only web interface"
}

# Parse command line arguments
case "${1:-start}" in
    "start"|"")
        start_services
        ;;
    "api-only")
        echo -e "${BLUE}Starting REST API only on http://localhost:8000${NC}"
        python3 api.py
        ;;
    "web-only")
        echo -e "${BLUE}Starting Web Interface only on http://localhost:8501${NC}"
        streamlit run app.py --server.port 8501
        ;;
    "help"|"-h"|"--help")
        show_help
        ;;
    *)
        echo -e "${RED}❌ Unknown option: $1${NC}"
        show_help
        exit 1
        ;;
esac