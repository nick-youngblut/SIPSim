#!/usr/bin/env python

# import
## batteries
import sys,os
import functools
import tempfile
## 3rd party
import numpy as np
## application libraries
libDir = os.path.dirname(__file__)
#libDir = os.path.join(scriptDir, '../lib/')
#sys.path.append(libDir)
import Utils



def parse_BD_windows(x):
    BDs = [map(float, y.split('-')) for y in x.split(',')]
    return BDs


def _HRSIP_by_window(BD_min, BD_max, prefix, Uargs):
    """HRSIP pipeline wrapper
    """
    # trimming physeq to just fractions in BD window
    exe = os.path.join(libDir, 'R', 'phyloseq_edit.r')
    inFile = Uargs['<phyloseq>']
    editFile = prefix + '_{}-{}.physeq'.format(BD_min, BD_max)
    editFile = os.path.split(editFile)[1]
    cmd = '{} --BD_min {} --BD_max {} {} > {}'
    cmd = cmd.format(exe, BD_min, BD_max,inFile, editFile)
    Utils.sys_call(cmd)

    # calling DESeq2                                                      
    exe = os.path.join(libDir, 'R', 'phyloseq_DESeq2.r')
    outFile = os.path.splitext(editFile)[0] + '_DESeq2'
    cmd = '{} --log2 {} --hypo {} --cont {} --treat {} --label {} --occur {}' \
          ' --padj {} {} > {}'
    cmd = cmd.format(exe, Uargs['--log2'], Uargs['--hypo'], Uargs['--cont'],
                     Uargs['--treat'], BD_min, Uargs['--occur'], 
                     Uargs['--padj'], editFile, outFile)
    Utils.sys_call(cmd)
    
    # returning output file
    return outFile


def HRSIP_by_window(BD_windows, prefix, Uargs):
    sys.stderr.write('# Running multi-window HR-SIP\n')
    
    res_files = []
    for BD_min,BD_max in BD_windows:        
        msg = '# HR-SIP on BD window: {}-{}\n'.format(BD_min, BD_max)
        sys.stderr.write(msg)
        resFile = _HRSIP_by_window(BD_min, BD_max, prefix, Uargs)
        res_files.append([resFile, BD_min, BD_max])
    return res_files


def _HRSIP_hierarchical(BD_min, BD_max, prefix, Uargs, tmpFile, i):
    """HRSIP pipeline wrapper
    """
    # trimming physeq to just fractions in BD window (and non-incorps)
    exe = os.path.join(libDir, 'R', 'phyloseq_edit.r')
    inFile = Uargs['<phyloseq>']
    editFile = prefix + '_{}-{}.physeq'.format(BD_min, BD_max)
    editFile = os.path.split(editFile)[1]
    if i < 1:
        cmd = '{} --BD_min {} --BD_max {} {} > {}'
        cmd = cmd.format(exe, BD_min, BD_max, inFile, editFile)
    else:
        cmd = '{} --BD_min {} --BD_max {} --taxa {} {} > {}'
        cmd = cmd.format(exe, BD_min, BD_max, tmpFile, inFile, editFile)
    Utils.sys_call(cmd)

    # calling DESeq2                                                      
    exe = os.path.join(libDir, 'R', 'phyloseq_DESeq2.r')
    outFile = os.path.splitext(editFile)[0] + '_DESeq2'
    cmd = '{} --log2 {} --hypo {} --cont {} --treat {} --label {}' \
          ' --occur {} --padj {} {} > {}'
    cmd = cmd.format(exe, Uargs['--log2'], Uargs['--hypo'], Uargs['--cont'],
                     Uargs['--treat'], BD_min, Uargs['--occur'], 
                     Uargs['--padj'], editFile, outFile)
    Utils.sys_call(cmd)

    # updating list of non-incorporators
    exe = os.path.join(libDir, 'R', 'DESeq2_listTaxa.r')
    cmd = '{} --inv --padj {} --log2 {} {} > {}'
    cmd = cmd.format(exe, Uargs['--padj'], Uargs['--log2'], outFile, tmpFile)
    Utils.sys_call(cmd)
    
    # returning output file
    return outFile


def HRSIP_hierarchical(BD_windows, prefix, Uargs):
    sys.stderr.write('# Running hierarchical multi-window HR-SIP\n')

    res_files = []
    tf = tempfile.NamedTemporaryFile()
    tmpFile = tf.name
    for i,(BD_min,BD_max) in enumerate(BD_windows):        
        msg = '# HR-SIP on BD window: {}-{}\n'.format(BD_min, BD_max)
        sys.stderr.write(msg)
        resFile = _HRSIP_hierarchical(BD_min, BD_max, prefix, Uargs, tmpFile, i)
        res_files.append([resFile, BD_min, BD_max])
    return res_files


def combine_DESeq(files, prefix):
    """Combining multiple DESeq2 objects; global p-value adjustment
    """
    exe = os.path.join(libDir, 'R', 'DESeq2_combine.r')
    outFile = prefix + '_DESeq2'
    outFile = os.path.split(outFile)[1]
    cmd = '{} {} > {}'.format(exe, ' '.join(files), outFile)
    Utils.sys_call(cmd)
    return outFile


def HR_SIP(Uargs):
    # file prefix
    if Uargs['--prefix'].lower() == 'none':
        prefix = os.path.splitext(Uargs['<phyloseq>'])[0]
    else:
        prefix = Uargs['--prefix']

    # parsing BD-windows
    BD_windows = parse_BD_windows(Uargs['-w'])

    # ordering by heaviest BD window (based on mean BD of range)
    if Uargs['--hier'] == True:
        BD_windows.sort(key=lambda x: -np.mean(x))
    
    # HR-SIP on each BD window
    if Uargs['--hier'] == True:
        res_files = HRSIP_hierarchical(BD_windows, prefix, Uargs)    
    else:
        res_files = HRSIP_by_window(BD_windows, prefix, Uargs)

    # combining DESeq2 files (if multiple)
    if len(res_files) > 1:
        outFile = combine_DESeq([x[0] for x in res_files], prefix)
    else:
        print 'Combining DESeq objects'
        outFile = res_files[0][0]
        outFileNew = prefix + '_DESeq2'
        outFileNew = os.path.split(outFileNew)[1]
        os.rename(outFile, outFileNew)

    # status
    print '\nFile written: {}'.format(outFile)

    
