.PHONY: all preprocess mapping timecourse spatial publish

R ?= Rscript

all: preprocess mapping timecourse spatial publish

preprocess:
	$(R) scripts/00_run_all.R --steps=preprocess

mapping:
	$(R) scripts/00_run_all.R --steps=mapping

timecourse:
	$(R) scripts/00_run_all.R --steps=timecourse

spatial:
	$(R) scripts/00_run_all.R --steps=spatial

publish:
	$(R) scripts/00_run_all.R --steps=publish
