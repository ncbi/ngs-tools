import json
import sys

# fields excluded from comparison
FIELDS_EXCLUDED = (
    "CMP_LINKAGE_GROUP", 
    "LINKAGE_GROUP", 
    "RAW_NAME", 
    'PLATFORM', 
    'CS_KEY', 
    'LANE', 
    'POSITION', 
    'X', 
    'Y', 
    'SPOT_IDS_FOUND', 
    'TILE',
    'CLIP_QUALITY_LEFT',
    'CLIP_QUALITY_RIGHT')

# fields corrected: 0-length read removed 
FIELDS_CORRECTED = (
    'READ_LEN', 
    'RD_FILTER', 
    'READ_FILTER', 
    'READ_TYPE', 
    'READ_SEG', 
    'READ_START')


def obj_cleanup(obj):
    for attr in FIELDS_EXCLUDED:
        if attr in obj:
            del obj[attr]
    read_len = obj['READ_LEN']
    if (isinstance(read_len, list)):
        index = -1
        for i, val in enumerate(read_len):
            if val == 0:
                index = i
        if index >= 0:                
            for attr in FIELDS_CORRECTED:
                f = obj[attr]
                if isinstance(f, list):
                    f.pop(index)
                    if (len(f) == 1 and attr != 'READ_SEG'):
                        obj[attr] = f[0]

def cmp_read(read1, read2):
    status = 0
    for f in read1:
        if f in read2:
            if (read1[f] != read2[f]):
                print(f + ": " + str(read1[f]) + " != " + str(read2[f]));   
                status = 1
        else:
            print(f + " does not exist");
            status = 1
    return status

# Opening JSON file
f1 = open(sys.argv[1])
data1 = json.load(f1)
f1.close()

f2 = open(sys.argv[2])
data2 = json.load(f2)
f2.close()
if len(data1) != len(data2):
    raise Exception("Different number of objects")
for read in data1:
    obj_cleanup(read)
for read in data2:
    obj_cleanup(read)

status = 0
for i, read in enumerate(data1):
    if cmp_read(read, data2[i]) != 0:
        status = 1

for i, read in enumerate(data2):
    if cmp_read(read, data1[i]) != 0:
        status = 1

sys.exit(status)
 

