#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2019 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

#
# Makefile for legacy results analysis
#

# global targets:

all: plotting

clean: clean_plotting


# plotting targets:

plotting:
	@python scripts/benchmark_plotter.py
	@python scripts/speedup_plotter.py
	@python scripts/speed_plotter.py

clean_plotting:
	@rm -rf benchmark_plots/pdf/*.pdf
	@rm -rf benchmark_plots/png/*.png
