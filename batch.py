import os
# for nn in range(500, 10000, 500):
for ss in [40, 400, 4000, 40000, 400000]:
	os.system("cargo run --release -- -i {0}".format(ss))