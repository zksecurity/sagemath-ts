#!/bin/bash
#
# Property Testing Runner: SageMath and PARI/GP vs sagemath-ts/parigp-ts
#
# This script runs the same operations in SageMath, PARI/GP, and TypeScript,
# then compares the results to verify they match exactly.
#
# Usage:
#   ./run-comparison.sh              # Run all tests
#   ./run-comparison.sh --generate   # Only generate reference output
#   ./run-comparison.sh --compare    # Only run TypeScript tests
#   ./run-comparison.sh --sage       # Run only SageMath tests
#   ./run-comparison.sh --pari       # Run only PARI/GP tests
#
# Requirements:
#   - SageMath (sage command) must be installed for arithmetic tests
#   - PARI/GP (gp command) must be installed for elliptic curve tests
#   - bun must be installed
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Transcript directories
PARI_TRANSCRIPTS="$SCRIPT_DIR/transcripts/pari"
SAGE_TRANSCRIPTS="$SCRIPT_DIR/transcripts/sage"
TS_TRANSCRIPTS="$SCRIPT_DIR/transcripts/typescript"

# Create transcript directories
mkdir -p "$PARI_TRANSCRIPTS"
mkdir -p "$SAGE_TRANSCRIPTS"
mkdir -p "$TS_TRANSCRIPTS"

# Parse arguments
GENERATE_ONLY=false
COMPARE_ONLY=false
SAGE_ONLY=false
PARI_ONLY=false
VERBOSE=false

for arg in "$@"; do
  case $arg in
    --generate)
      GENERATE_ONLY=true
      shift
      ;;
    --compare)
      COMPARE_ONLY=true
      shift
      ;;
    --sage)
      SAGE_ONLY=true
      shift
      ;;
    --pari)
      PARI_ONLY=true
      shift
      ;;
    --verbose|-v)
      VERBOSE=true
      shift
      ;;
    *)
      ;;
  esac
done

# Check if SageMath is installed
check_sage() {
  # Try common sage paths
  SAGE_PATH=""
  for path in /usr/local/bin/sage /opt/homebrew/bin/sage /usr/bin/sage sage; do
    if command -v "$path" &> /dev/null; then
      SAGE_PATH="$path"
      break
    fi
  done

  if [ -z "$SAGE_PATH" ]; then
    echo -e "${YELLOW}Warning: SageMath (sage command) not found.${NC}"
    echo ""
    echo "To install SageMath:"
    echo ""
    echo "  macOS (Homebrew):"
    echo "    brew install --cask sage"
    echo ""
    echo "  Ubuntu/Debian:"
    echo "    sudo apt-get install sagemath"
    echo ""
    echo "  Or download from: https://www.sagemath.org/download.html"
    echo ""
    echo "After installation, run this script again."
    return 1
  fi
  echo "$SAGE_PATH"
}

# Check if PARI/GP is installed
check_pari() {
  # Use command -v to bypass aliases
  if ! command -v /usr/local/bin/gp &> /dev/null && ! command -v /opt/homebrew/bin/gp &> /dev/null; then
    # Try to find gp in common locations
    GP_PATH=""
    for path in /usr/local/bin/gp /opt/homebrew/bin/gp /usr/bin/gp; do
      if [ -x "$path" ]; then
        GP_PATH="$path"
        break
      fi
    done

    if [ -z "$GP_PATH" ]; then
      echo -e "${YELLOW}Warning: PARI/GP (gp command) not found.${NC}"
      echo ""
      echo "To install PARI/GP:"
      echo ""
      echo "  macOS (Homebrew):"
      echo "    brew install pari"
      echo ""
      echo "  Ubuntu/Debian:"
      echo "    sudo apt-get install pari-gp"
      echo ""
      echo "  Fedora/RHEL:"
      echo "    sudo dnf install pari-gp"
      echo ""
      echo "  Windows (WSL):"
      echo "    sudo apt-get install pari-gp"
      echo ""
      echo "After installation, run this script again."
      return 1
    fi
    echo "$GP_PATH"
  else
    # Return the path to gp
    for path in /opt/homebrew/bin/gp /usr/local/bin/gp; do
      if [ -x "$path" ]; then
        echo "$path"
        return 0
      fi
    done
  fi
}

# Run SageMath scripts
run_sage_tests() {
  echo -e "${BLUE}=== Running SageMath Tests ===${NC}"

  SAGE_CMD=$(check_sage)
  if [ $? -ne 0 ]; then
    echo "$SAGE_CMD"  # This contains the error message
    return 1
  fi

  echo -e "Using SageMath: ${GREEN}$SAGE_CMD${NC}"
  echo ""

  # Run arithmetic tests
  if [ -f "$SCRIPT_DIR/sage/arith.sage" ]; then
    echo -e "Running: ${YELLOW}arith.sage${NC}"
    "$SAGE_CMD" "$SCRIPT_DIR/sage/arith.sage" > "$SAGE_TRANSCRIPTS/arith.txt" 2>&1

    if [ $? -eq 0 ]; then
      echo -e "  ${GREEN}OK${NC} - Output saved to transcripts/sage/arith.txt"
      # Count results
      local result_count=$(grep -c "^RESULT:" "$SAGE_TRANSCRIPTS/arith.txt" 2>/dev/null || echo "0")
      echo -e "  Found ${CYAN}$result_count${NC} test results"
    else
      echo -e "  ${RED}FAILED${NC} - Check $SAGE_TRANSCRIPTS/arith.txt for errors"
      cat "$SAGE_TRANSCRIPTS/arith.txt"
      return 1
    fi
  else
    echo -e "  ${YELLOW}SKIP${NC} - arith.sage not found"
  fi

  echo ""
  return 0
}

# Run PARI/GP scripts
run_pari_tests() {
  echo -e "${BLUE}=== Running PARI/GP Tests ===${NC}"

  GP_CMD=$(check_pari)
  if [ $? -ne 0 ]; then
    echo "$GP_CMD"  # This contains the error message
    return 1
  fi

  echo -e "Using PARI/GP: ${GREEN}$GP_CMD${NC}"
  echo ""

  # Run elliptic curve tests
  if [ -f "$SCRIPT_DIR/pari/elliptic.gp" ]; then
    echo -e "Running: ${YELLOW}elliptic.gp${NC}"
    "$GP_CMD" -q < "$SCRIPT_DIR/pari/elliptic.gp" > "$PARI_TRANSCRIPTS/elliptic.txt" 2>&1

    if [ $? -eq 0 ]; then
      echo -e "  ${GREEN}OK${NC} - Output saved to transcripts/pari/elliptic.txt"
      # Count results
      local result_count=$(grep -c "^RESULT:" "$PARI_TRANSCRIPTS/elliptic.txt" 2>/dev/null || echo "0")
      echo -e "  Found ${CYAN}$result_count${NC} test results"
    else
      echo -e "  ${RED}FAILED${NC} - Check $PARI_TRANSCRIPTS/elliptic.txt for errors"
      return 1
    fi
  else
    echo -e "  ${YELLOW}SKIP${NC} - elliptic.gp not found"
  fi

  echo ""
  return 0
}

# Run TypeScript tests for SageMath comparison
run_sage_typescript_tests() {
  echo -e "${BLUE}=== Running TypeScript Tests (SageMath Comparison) ===${NC}"

  cd "$PROJECT_ROOT"

  # Check if SageMath output exists
  if [ ! -f "$SAGE_TRANSCRIPTS/arith.txt" ]; then
    echo -e "${YELLOW}Warning: SageMath output not found. Run with --generate or --sage first.${NC}"
    echo ""
  fi

  # Run the bun tests for arithmetic comparison
  if [ -f "$SCRIPT_DIR/typescript/arith-comparison.test.ts" ]; then
    echo -e "Running: ${YELLOW}arith-comparison.test.ts${NC}"

    if bun test "$SCRIPT_DIR/typescript/arith-comparison.test.ts" 2>&1; then
      echo -e "  ${GREEN}OK${NC} - All arithmetic tests passed"
    else
      echo -e "  ${RED}FAILED${NC} - Some arithmetic tests failed"
      return 1
    fi
  else
    echo -e "  ${YELLOW}SKIP${NC} - arith-comparison.test.ts not found"
  fi

  echo ""
  return 0
}

# Run TypeScript tests for PARI/GP comparison
run_pari_typescript_tests() {
  echo -e "${BLUE}=== Running TypeScript Tests (PARI/GP Comparison) ===${NC}"

  cd "$PROJECT_ROOT"

  # Check if PARI output exists
  if [ ! -f "$PARI_TRANSCRIPTS/elliptic.txt" ]; then
    echo -e "${YELLOW}Warning: PARI output not found. Run with --generate or --pari first.${NC}"
    echo ""
  fi

  # Run the bun tests
  if [ -f "$SCRIPT_DIR/typescript/elliptic.test.ts" ]; then
    echo -e "Running: ${YELLOW}elliptic.test.ts${NC}"

    if bun test "$SCRIPT_DIR/typescript/elliptic.test.ts" 2>&1; then
      echo -e "  ${GREEN}OK${NC} - All elliptic curve tests passed"
    else
      echo -e "  ${RED}FAILED${NC} - Some elliptic curve tests failed"
      return 1
    fi
  else
    echo -e "  ${YELLOW}SKIP${NC} - elliptic.test.ts not found"
  fi

  echo ""
  return 0
}

# Compare transcripts and show summary
compare_results() {
  echo -e "${BLUE}=== Results Summary ===${NC}"

  local all_passed=true

  # Check SageMath results
  if [ -f "$SAGE_TRANSCRIPTS/arith.txt" ]; then
    # Extract just the RESULT lines for comparison
    grep "^RESULT:" "$SAGE_TRANSCRIPTS/arith.txt" > "$SAGE_TRANSCRIPTS/arith.results.txt" 2>/dev/null || true

    echo -e "SageMath arithmetic results: ${YELLOW}transcripts/sage/arith.results.txt${NC}"
    local sage_count=$(wc -l < "$SAGE_TRANSCRIPTS/arith.results.txt" | tr -d ' ')
    echo -e "  Found ${GREEN}$sage_count${NC} test results"
  else
    echo -e "${YELLOW}No SageMath output available${NC}"
  fi

  # Check PARI/GP results
  if [ -f "$PARI_TRANSCRIPTS/elliptic.txt" ]; then
    # Extract just the RESULT lines for comparison
    grep "^RESULT:" "$PARI_TRANSCRIPTS/elliptic.txt" > "$PARI_TRANSCRIPTS/elliptic.results.txt" 2>/dev/null || true

    echo -e "PARI/GP elliptic results: ${YELLOW}transcripts/pari/elliptic.results.txt${NC}"
    local pari_count=$(wc -l < "$PARI_TRANSCRIPTS/elliptic.results.txt" | tr -d ' ')
    echo -e "  Found ${GREEN}$pari_count${NC} test results"
  else
    echo -e "${YELLOW}No PARI/GP output available${NC}"
  fi

  echo ""
  return 0
}

# Show help
show_help() {
  echo "Property Testing Runner: SageMath/PARI vs sagemath-ts/parigp-ts"
  echo ""
  echo "Usage:"
  echo "  $0              Run all tests (generate output + run TS tests)"
  echo "  $0 --generate   Only generate reference output (SageMath + PARI)"
  echo "  $0 --compare    Only run TypeScript comparison tests"
  echo "  $0 --sage       Run only SageMath arithmetic tests"
  echo "  $0 --pari       Run only PARI/GP elliptic curve tests"
  echo "  $0 --verbose    Show detailed output"
  echo "  $0 --help       Show this help"
  echo ""
  echo "Requirements:"
  echo "  - SageMath (sage command) - for arithmetic function tests"
  echo "  - PARI/GP (gp command) - for elliptic curve tests"
  echo "  - bun - for running TypeScript tests"
  echo ""
  echo "Test modules:"
  echo "  - sage/arith.sage         -> Arithmetic functions (gcd, factor, etc.)"
  echo "  - pari/elliptic.gp        -> Elliptic curve operations"
  echo ""
  echo "Elliptic curve test categories (pari/elliptic.gp):"
  echo "  - ellinit: Curve initialization (3 tests)"
  echo "  - ellcard: Cardinality computation (8 tests)"
  echo "  - Point operations: elladd, ellmul, ellneg (7 tests)"
  echo "  - ellorder: Point order (3 tests)"
  echo "  - Identity element tests (5 tests)"
  echo "  - Discriminant and j-invariant (2 tests)"
  echo "  - ellisoncurve: Curve membership (2 tests)"
  echo "  - ellgroup: Group structure (3 tests)"
  echo "  - Associativity tests (2 tests)"
  echo "  - elllog: Discrete logarithm (5+ tests)"
  echo "  - elldivpol: Division polynomials (7 tests)"
  echo "  - elltatepairing: Tate pairing (4 tests)"
  echo "  - ellweilpairing: Weil pairing (3 tests)"
  echo "  - Known curves: secp256k1-like, P-256-like (3 tests)"
  echo "  - Random curves (3 tests)"
  echo "  - Edge cases: Small primes, singular curves (4 tests)"
  echo "  - Large field tests (2 tests)"
  echo "  - Pairing bilinearity (2 tests)"
}

# Main
main() {
  if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    show_help
    exit 0
  fi

  echo ""
  echo -e "${BLUE}╔═══════════════════════════════════════════════════════════════════╗${NC}"
  echo -e "${BLUE}║   Property Testing: SageMath/PARI vs sagemath-ts/parigp-ts       ║${NC}"
  echo -e "${BLUE}╚═══════════════════════════════════════════════════════════════════╝${NC}"
  echo ""

  local exit_code=0

  if $SAGE_ONLY; then
    run_sage_tests || exit_code=1
    if [ $exit_code -eq 0 ]; then
      run_sage_typescript_tests || exit_code=1
    fi
  elif $PARI_ONLY; then
    run_pari_tests || exit_code=1
    if [ $exit_code -eq 0 ]; then
      run_pari_typescript_tests || exit_code=1
    fi
  elif $GENERATE_ONLY; then
    run_sage_tests || exit_code=1
    run_pari_tests || exit_code=1
  elif $COMPARE_ONLY; then
    run_sage_typescript_tests || exit_code=1
    run_pari_typescript_tests || exit_code=1
  else
    # Run everything
    run_sage_tests || exit_code=1
    run_pari_tests || exit_code=1
    if [ $exit_code -eq 0 ]; then
      run_sage_typescript_tests || exit_code=1
      run_pari_typescript_tests || exit_code=1
    fi
  fi

  compare_results

  if [ $exit_code -eq 0 ]; then
    echo -e "${GREEN}=== All tests completed successfully ===${NC}"
  else
    echo -e "${RED}=== Some tests failed ===${NC}"
  fi

  exit $exit_code
}

main "$@"
