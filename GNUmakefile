SOURCE_DIR := committor
TEST_DIR := tests
SCRIPTS_DIR := scripts

SOURCES := $(wildcard $(SOURCE_DIR)/*.py)
TESTS := $(wildcard $(TEST_DIR)/*.py)

AUTO_GENERATED_SOURCES := $(SOURCE_DIR)/version.py

all: $(AUTO_GENERATED_SOURCES)

$(SOURCE_DIR)/version.py:
	@scripts/make-version > $@

.PHONY: all $(SOURCE_DIR)/version.py
