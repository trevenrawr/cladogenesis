import os
# for nn in range(500, 10000, 500):
for ss in [1.0, 10.0, 100.0, 1000.0, 10000.0]:
	os.system("cargo run --release -- -m {0}".format(ss))