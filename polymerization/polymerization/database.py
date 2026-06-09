# functions to re-index giggle gene queries
import os
import subprocess
import sys
import gzip
import tempfile
def prep_giggle_clean(path_giggle_clean, outfile, sample_column=10, header=True):
    try:
        with gzip.open(path_giggle_clean, 'rt') as f_in, tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_remainder:
            header_lines = []
            for line in f_in:
                if line.startswith("#") and header:
                    header_lines.append(line)
                else:
                    temp_remainder.write(line)
            temp_remainder.close()

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_combined:
            temp_combined.writelines(header_lines)
            temp_combined.flush()
            subprocess.run(f"cut --complement -f {sample_column} {temp_remainder.name} >> {temp_combined.name}", shell=True, check=True)
            subprocess.run(f"bgzip -c {temp_combined.name} > {outfile}", shell=True, check=True)

        os.remove(temp_remainder.name)
        os.remove(temp_combined.name)
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False