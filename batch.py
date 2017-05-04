import os

for rr in [0.00005, 0.0005, 0.005]:
	os.system("cargo run --release -- -p {} -a false -b 100".format(rr))