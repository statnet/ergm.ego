data_objects=fmhfit

data_files=$(addprefix ../data/,$(data_objects:=.rda))

default: $(data_files)

$(data_files): ../data/%.rda: %.R
	Rscript $<

.PHONY: default

