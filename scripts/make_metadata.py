#!/usr/env python

import sys
import os
import shutil
import tarfile

tmpdir = sys.argv[1]
zipname = sys.argv[2]
files = sys.argv[3:]
samples = [f.split("/")[-3] for f in files]



print("name,timepoint,biorep,path")



for f, s in zip(files, samples):
    r = None

    dest = os.path.join("kallisto", s, "abundance.tsv")
    tmpfile = os.path.join(tmpdir, dest)
    
    if not os.path.exists(os.path.split(tmpfile)[0]):
        os.makedirs(str(os.path.split(tmpfile)[0]))
    shutil.copyfile(f,tmpfile)

    #now copy also needed files abundance.h5 and run_info.json
    h5i = f.replace(".tsv", ".h5")
    h5o = os.path.join(tmpdir, "kallisto", s, "abundance.h5")
    shutil.copyfile(h5i, h5o)

    jsi = f.replace("abundance.tsv", "run_info.json")
    jso = os.path.join(tmpdir, "kallisto", s, "run_info.json")
    shutil.copyfile(jsi,jso)

    md = s.split("_")
    r = "{},{},{},{}".format(md[0], md[1], md[2], dest)
    print(r)


if not os.path.exists("results"):
    os.mkdir("results")

src = os.path.join(os.sep, tmpdir, "kallisto")
with tarfile.open(zipname, "w:gz") as tar:
    tar.add(src,arcname=os.path.basename(src))
shutil.rmtree(src)