# The "make" targets are:
#    wheel: build a Python wheel in "dist" directory
#    install: build wheel and install
#    test: run ChimeraX with test commands
#    debug: run ChimeraX with debugging flag set
#    clean: remove files used in building wheel

# Theoretically, you should only need to change
# 1. the CHIMERAX_APP definition (or define in environment),
# 2. the "devel build" and "devel install" commands in the
#    "install" and "wheel" targets, (e.g., if you want to
#    install into ChimeraX.app instead of for the current
#    user, add "user false"), and
# 3. the "generated_files" and "clean_generated_files" targets
#    (e.g., if you need to generate files such as self-signed
#    certificates or data files from template expansion)

# Define where ChimeraX is installed.
OS=$(shell uname -s)
# We're on Windows
ifeq ($(filter $(OS),Linux Darwin),)
OS=$(shell uname -o)
ifneq ($(filter $(OS),Cygwin Msys),)
OS=Windows
endif
endif

# CHIMERAX_APP is the ChimeraX install folder
# We use ?= to let CHIMERAX_APP environment variable override
ifeq ($(OS),Windows)
# Windows
ifndef RELEASE
CHIMERAX_APP ?= "c:/Program Files/ChimeraX_Daily"
else
CHIMERAX_APP ?= "c:/Program Files/ChimeraX"
endif
endif

ifeq ($(OS),Darwin)
# Mac
ifndef RELEASE
CHIMERAX_APP ?= "/Users/greg/Desktop/ChimeraX-1.6.1.app/"
else
CHIMERAX_APP ?= "/Users/greg/Desktop/ChimeraX-1.6.1.app/"
endif
endif
# Platform-dependent settings.  Should not need fixing.
# For Windows, we assume Cygwin is being used.
ifeq ($(OS),Windows)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX.exe
endif
ifeq ($(OS),Darwin)
export MACOSX_DEPLOYMENT_TARGET=10.13
CHIMERAX_EXE = $(CHIMERAX_APP)/Contents/bin/ChimeraX
endif
ifeq ($(OS),Linux)
ifndef RELEASE
CHIMERAX_EXE = $(shell which chimerax-daily)
else
CHIMERAX_EXE = $(shell which chimerax)
endif
endif
RUN = $(CHIMERAX_EXE) --nogui --exit --cmd

PYSRCS = $(wildcard src/*.py)
CSRCS = $(wildcard src/*.cpp)
SRCS = $(PYSRCS) $(CSRCS)

# If you want to install into ChimeraX.app, add "user false"
# to the "devel build" and "devel install" commands.
# By default, we install for just the current user.

install:	pyproject.toml $(SRCS) generated_files
	$(RUN) "devel install . exit true"

wheel:	pyproject.toml $(SRCS) generated_files
	$(RUN) "devel build . exit true"

docs:
	$(CHIMERAX_EXE) -m sphinx docs/source src/docs/user

test:
	for t in $(wildcard test*.cxc) $(wildcard test*.py);\
		do $(CHIMERAX_EXE) --exit --nogui $$t;\
	done

debug:
	$(CHIMERAX_EXE) --debug

clean:	clean_generated_files
	if [ -x $(CHIMERAX_EXE) ]; then \
		$(RUN) "devel clean . exit true" ; \
	else \
		rm -rf build dist *.egg-info src/__pycache__ ; \
	fi

pylint:
	$(CHIMERAX_EXE) -m flake8 $(filter %.py, $(SRCS))

# Modify targets below if you have files that need
# to be generated before building the bundle and/or
# removed when cleaning the bundle

generated_files:
	# Generate files

clean_generated_files:
	# Remove generated files

.PHONY: docs
