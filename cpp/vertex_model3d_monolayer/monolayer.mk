# monolayer.mk
SHELL := /bin/bash

# Base include directories for the monolayer project
INCS := $(shell \
  base="-I./src/ -I./examples/ -I./build/ -I./src/eigen-3.4.0"; \
  for dep in boost@1.85 libomp gsl; do \
    prefix=`brew --prefix $$dep`; \
    base="$$base -I$$prefix/include"; \
  done; \
  echo $$base \
)

# Library paths for the monolayer project
LIBPATH := $(shell \
  path=""; \
  for dep in boost@1.85 libomp gsl; do \
    prefix=`brew --prefix $$dep`; \
    path="$$path -L$$prefix/lib"; \
  done; \
  echo $$path \
)
