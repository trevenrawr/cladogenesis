import os
for _ in range(1000):
	os.system("cargo run --release -- -r false -a false".format(ss))