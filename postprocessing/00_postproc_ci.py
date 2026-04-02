import sys
sys.path.append("/glade/u/home/tomasc/repos/pynanigans") # Add pynanigans to PYTHONPATH (HPC path; pip-installed in CI)

print("Starting CI post-processing")

#+++ Define run options
simdata_path = "../simulations/data/"
simname_base = "ci_test"
runs = [{}]  # single run, no parameter suffix
#---

exec(open("01_create_aaaa.py").read())
exec(open("02_create_xyza.py").read())
