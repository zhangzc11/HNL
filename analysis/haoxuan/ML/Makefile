# Makefile for HNL KS plotting workflow
# Usage examples:
#  make help
#  make extract           # build filtered_branches.csv
#  make ks                # run ks_branch_comparison.py
#  make plot_mass MASS=5GeV   # plot for fixed mass
#  make plot_tau TAU=0ps      # plot for fixed tau
#  make plot_both MASS=5GeV TAU=0ps
#  make show              # list generated pngs
#  make clean             # remove generated pngs

PY ?= python
ROOT_PY ?= /home/haoxuan/anaconda3/envs/root_env/bin/python
CSV_DIR = ./csv
TOP30 = $(CSV_DIR)/top30_ks_branches.csv
FILTERED = $(CSV_DIR)/filtered_branches.csv
MC_DIR = /pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/
DATA_FILE = data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root
PLOT_SCRIPT = plot_ks_images.py
KS_SCRIPT = ks_branch_comparison.py
EXTRACT_SCRIPT = extract_common_branches.py
PNG_DIR = ./png

.PHONY: help extract ks plot_mass plot_tau plot_both show clean env_check plot_all_tau plot_all_mass plot_all_both

help:
	@echo "Makefile targets for HNL plotting"
	@echo "  make extract            -> run $(EXTRACT_SCRIPT) to produce $(FILTERED)"
	@echo "  make ks                 -> run $(KS_SCRIPT) to produce $(TOP30)"
	@echo "  make plot_mass MASS=5GeV-> run $(PLOT_SCRIPT) fixed mass (requires ROOT env)"
	@echo "  make plot_tau TAU=0ps   -> run $(PLOT_SCRIPT) fixed tau (requires ROOT env)"
	@echo "  make plot_both MASS=5GeV TAU=0ps -> run both plots"
	@echo "  make show               -> list generated png files"
	@echo "  make clean              -> remove generated png files"
	@echo "Note: set ROOT_PY to your ROOT-enabled python if different."

# Run branch extraction (uses uproot)
extract: $(FILTERED)

$(FILTERED): $(EXTRACT_SCRIPT)
	@echo "Running branch extraction..."
	$(PY) $(EXTRACT_SCRIPT)
	@echo "-> $(FILTERED) should be created (or check script output)"

# Run KS comparison (requires filtered branches and data/MC access)
ks: $(TOP30)

$(TOP30): $(KS_SCRIPT) $(FILTERED)
	@echo "Running KS comparison..."
	$(PY) $(KS_SCRIPT)
	@echo "-> $(TOP30) generated"

# Plot for fixed mass (uses ROOT python)
plot_mass: $(TOP30)
	@if [ -z "$(MASS)" ]; then echo "Please provide MASS, e.g. make plot_mass MASS=5GeV"; exit 1; fi
	@echo "Plotting for Mass=$(MASS) using ROOT python ($(ROOT_PY))..."
	$(ROOT_PY) $(PLOT_SCRIPT) 1 $(MASS)
	@echo "Plots saved under $(PNG_DIR)/M$(MASS)/"

plot_all_mass: $(TOP30)
	@echo "Plotting for all masses using ROOT python ($(ROOT_PY))..."
	$(ROOT_PY) $(PLOT_SCRIPT) 1 5GeV
	$(ROOT_PY) $(PLOT_SCRIPT) 1 10GeV
	$(ROOT_PY) $(PLOT_SCRIPT) 1 15GeV
	$(ROOT_PY) $(PLOT_SCRIPT) 1 20GeV
	$(ROOT_PY) $(PLOT_SCRIPT) 1 30GeV
	$(ROOT_PY) $(PLOT_SCRIPT) 1 50GeV
	@echo "Plots saved under $(PNG_DIR)/all_mass/"

# Plot for fixed tau (uses ROOT python)
plot_tau: $(TOP30)
	@if [ -z "$(TAU)" ]; then echo "Please provide TAU, e.g. make plot_tau TAU=0ps"; exit 1; fi
	@echo "Plotting for Tau=$(TAU) using ROOT python ($(ROOT_PY))..."
	$(ROOT_PY) $(PLOT_SCRIPT) 2 $(TAU)
	@echo "Plots saved under $(PNG_DIR)/tau$(TAU)/"

plot_all_tau: $(TOP30)
	@echo "Plotting for all taus using ROOT python ($(ROOT_PY))..."
	$(ROOT_PY) $(PLOT_SCRIPT) 2 0ps
	$(ROOT_PY) $(PLOT_SCRIPT) 2 10ps
	$(ROOT_PY) $(PLOT_SCRIPT) 2 100ps
	$(ROOT_PY) $(PLOT_SCRIPT) 2 1000ps
	@echo "Plots saved under $(PNG_DIR)/all_tau/"

plot_both: plot_mass plot_tau
plot_all_both: plot_all_mass plot_all_tau

# List generated pngs
show:
	@echo "PNG directories:"
	@ls -1 $(PNG_DIR) || echo "No png directory"
	@echo "Contents:"
	@ls -R $(PNG_DIR) || true

# Clean generated pngs
clean:
	@echo "Removing png directory $(PNG_DIR)"
	rm -rf $(PNG_DIR)
	@echo "Done."

# Quick check for ROOT python
env_check:
	@echo "Checking ROOT python at $(ROOT_PY)"
	@$(ROOT_PY) -c "import sys; import importlib\ntry:\n    importlib.import_module('ROOT'); print('ROOT available')\nexcept Exception as e:\n    print('ROOT import failed:', e); sys.exit(2)"

# Default target
all: help
