import os
import sys
import glob
import re
import argparse
import logging
import pandas as pd
import contextlib # <--- NEW IMPORT
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# Import your pipeline logic
from batch_driver import run_sightline_pipeline, SPECIES_CONFIG

def scan_for_ids(data_dir, species):
    # ... (Same as before) ...
    if species not in SPECIES_CONFIG:
        raise ValueError(f"Unknown species: {species}")
    
    pattern = f"{species}_doublet_los_*_snr10.txt"
    search_path = os.path.join(data_dir, pattern)
    files = glob.glob(search_path)
    ids = []
    regex = re.compile(r"los_(\d+)_snr10")
    for f in files:
        match = regex.search(f)
        if match: ids.append(match.group(1))
    return sorted(list(set(ids)))

def worker_wrapper(args):
    """
    Wrapper to unpack arguments and NUKE all output.
    """
    los_id, data_dir, out_dir, species, trim = args
    
    # 1. Nuclear Option: Disable all logging
    # This prevents ANY log message (CRITICAL and below) from being processed.
    # Since this runs in a separate process, it won't affect the main orchestrator.
    logging.disable(logging.CRITICAL)

    # 2. Redirect stdout/stderr
    # This catches "print()" calls and C-level writes (if supported by OS/buffer)
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            try:
                # Run Pipeline
                result_path = run_sightline_pipeline(data_dir, out_dir, los_id, species, trim)
                
                if result_path and os.path.exists(result_path):
                    return {"id": los_id, "status": "success", "path": result_path}
                else:
                    return {"id": los_id, "status": "failed", "error": "No output generated"}
                    
            except Exception as e:
                # We return the error message string so the Main Orchestrator 
                # can print it nicely in the summary/progress bar.
                return {"id": los_id, "status": "error", "error": str(e)}

def aggregate_results(output_dir, species, final_catalog_name="master_catalog.csv"):
    # ... (Same as before) ...
    print(f"\nAggregating results for {species}...")
    all_files = glob.glob(os.path.join(output_dir, f"results_los_*_{species}.csv"))
    
    if not all_files:
        print("No result files found to aggregate.")
        return

    combined_df = []
    for f in tqdm(all_files, desc="Merging CSVs"):
        try:
            df = pd.read_csv(f, comment='#')
            basename = os.path.basename(f)
            match = re.search(r"los_(\d+)_", basename)
            if match: df['los_id'] = match.group(1)
            combined_df.append(df)
        except Exception as e:
            pass # Silent skip for aggregation

    if combined_df:
        master = pd.concat(combined_df, ignore_index=True)
        cols = ['los_id'] + [c for c in master.columns if c != 'los_id']
        master = master[cols]
        save_path = os.path.join(output_dir, final_catalog_name)
        master.to_csv(save_path, index=False)
        print(f"Master catalog saved to: {save_path}")
        return save_path
    return None

def main():
    parser = argparse.ArgumentParser(description="Parallel Orchestrator")
    parser.add_argument('--data', type=str, required=True)
    parser.add_argument('--out', type=str, default="batch_output")
    parser.add_argument('--species', type=str, default='OVI')
    parser.add_argument('--cores', type=int, default=4)
    parser.add_argument('--trim', type=float, default=500.0)
    args = parser.parse_args()
    
    print(f"Scanning {args.data} for {args.species}...")
    los_ids = scan_for_ids(args.data, args.species)
    print(f"Found {len(los_ids)} lines of sight.")
    if not los_ids: sys.exit("No input files found.")

    os.makedirs(args.out, exist_ok=True)

    tasks = [(uid, args.data, args.out, args.species, args.trim) for uid in los_ids]
    results = []
    
    print(f"Starting processing on {args.cores} cores...")
    
    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = {executor.submit(worker_wrapper, task): task[0] for task in tasks}
        
        # Now the progress bar will be clean!
        with tqdm(total=len(los_ids), desc="Processing LOS") as pbar:
            for future in as_completed(futures):
                res = future.result()
                results.append(res)
                
                if res['status'] != 'success':
                    # We manually print errors *outside* the silenced context
                    pbar.write(f"LOS {res['id']} failed: {res.get('error')}")
                
                pbar.update(1)

    success_count = sum(1 for r in results if r['status'] == 'success')
    print(f"\nProcessing Complete. Success: {success_count}/{len(los_ids)}")
    aggregate_results(args.out, args.species)

if __name__ == "__main__":
    os.environ['MPLBACKEND'] = 'Agg'
    main()