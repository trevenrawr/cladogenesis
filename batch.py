import os
# for nn in range(500, 10000, 500):
for ss in [1000, 2000, 5000, 8000, 13000]:
	os.system("cargo run --release -- -n {0} -a false".format(ss))