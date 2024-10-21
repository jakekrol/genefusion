import os,sys
import time
from multiprocessing import Pool
def stix(index, db, chr, left_start, left_end, right_start, right_end, outfile,type,chdir=None):
    if chdir:
        os.chdir(chdir)
    cmd = f"stix -i {index} -d {db} -t {type} -l {chr}:{left_start}-{left_end} -r {chr}:{right_start}-{right_end} -s 500 > {outfile}"
    print(cmd)
    os.system(cmd)
def stix_sharded(dir_shard, chr, left_start, left_end, right_start, right_end,outdir,outfile,type,processes,index_name='index',db_name='stix.ped.db'):
    print('begin sharded stix')
    processes = int(processes)
    args = locals()
    for key, value in args.items():
        print(f"{key}: {value}")
    time.sleep(3)
    shards = os.listdir(dir_shard)
    if not shards:
        print('no shards found at',dir_shard)
        return
    shards = [os.path.join(dir_shard,shard) for shard in shards]
    data = []
    for shard in shards:
         shard_out = os.path.join(outdir,f'{os.path.basename(shard)}.{outfile}')
         data.append((index_name,db_name,chr,left_start,left_end,right_start,right_end,shard_out,type,shard))
    print('data','\n',data[0])
    print('performing stix queries over shards at',dir_shard, 'with',processes,'processes')
    with Pool(processes=processes) as pool:
        pool.starmap(stix, data)
    print('end sharded stix')
