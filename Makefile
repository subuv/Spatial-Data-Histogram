MODULES = spatialDataHistogram
PGXS := $(shell pg_config --pgxs)
include $(PGXS)