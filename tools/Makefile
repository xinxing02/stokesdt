include ../makevars.in

SUBDIR = calc_diff mol_gen

all:
	@for dir in $(SUBDIR); \
          do \
            make -C $${dir} || exit 1; \
          done

install:
	@for dir in $(SUBDIR); \
          do \
            make -C $${dir} install || exit 1; \
          done

clean:
	@for dir in $(SUBDIR); \
          do \
            make -C $${dir} clean || exit 1; \
          done
