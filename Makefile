# Usage:
# make	# compile
# make clean  # remove object and executable files

.PHONY = all


all:
	$(MAKE) -C thrid_party/libMeshb/
	$(MAKE) -C src/
	

clean:
	$(MAKE) -C thrid_party/libMeshb/ 'clean'
	$(MAKE) -C src/ 'clean'