the_headers = oratio.h \
	clapack.h error.h \
	matrix.h \
	vector.h \
	bsplines.h \
	cblas.h \
	cblaswr_cxml.h \
	constants.h \
	ctools.h ctools_error.h ctools_macros.h cfile.h cstring.h \
	csystem.h ctime.h ctimer.h \
	getopt.h \
	exception.h \
	dsarray.h dsbtree.h \
	dsconvert.h \
	dsqueue.h \
	dstruct.h ds_functions.h ds_sarray.h ds_search.h ds_sort.h \
	dslist.h dsslist.h \
	empirical.h \
	f2c.h \
	faure.h \
	functions.h \
	generation.h \
	histogram.h \
	infinity.h \
	log.h \
	lognormal.h \
	memory.h \
	mlhs.h \
	moca.h halton.h \
	moca_io.h \
	natural.h \
	norm.h \
	ntools.h nt_coordinates.h nt_derive.h \
	nt_gamma.h nt_io.h nt_polynomial.h \
	nt_simplex.h nt_utilities.h \
	percentile.h \
	ran.h ran_rng.h ran_bernouilli.h ran_copula.h \
	ran_mvn.h ran_walk.h \
	ran_brownian.h ran_mpar.h \
	sim.hpp  sim_accumulate.hpp \
	sobol.h \
	stat.h \
	st_descriptive.h st_likelihood.h st_obs.h st_utils.h \
	st_betasymmetrical.h st_cauchy.h st_chi2.h st_discrete.h\
	st_exponential.h \
	st_fisher.h st_geometric.h  st_gumbel.h st_laplace.h st_logistic.h \
	st_lognormal.h st_normal.h \
	st_pareto.h st_poisson.h \
	st_rayleigh.h st_student.h st_triangular.h st_uniform.h st_weibull.h \
	st_mvn.h \
	st_tally.h \
	stat_empirical.h \
	stat_error.h \
	symbol.h \
	truncation.h \
	ttest.h \
	values.h \
	graph.h

the_src = credits.c

lib_LTLIBRARIES = liboratio.la

liboratio_la_SOURCES = $(the_headers) $(the_src)

liboratio_la_LIBADD =
liboratio_la_LDFLAGS = -no-undefined

oratioincludedir = $(pkgincludedir)
oratioinclude_HEADERS = $(the_headers)

EXTRA_DIST = $(the_headers)
