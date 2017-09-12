#!/usr/bin/env python

import os, sys
import subprocess as sp
import glob
import argparse
import ROOT as rt

from SMSConfig import sms_models, VERSION

def do_command(cmd, no_exec=True):
    """
    Input: list of strings that should be joined to form the command.
    """
    print ' '.join(cmd)
    if not no_exec:
        sp.call(cmd)

def get_limit_dir(model):
    return '/eos/cms/store/group/phys_susy/razor/Run2Analysis/Limits/RazorInclusive2016/{}/{}'.format(
            VERSION, model)

def write_done_file(model, box):
    """
    Creates text file containing filenames of completed jobs.
    Returns the name of the file.
    """
    done_jobs = glob.glob(get_limit_dir(model)
            +'/higgsCombineMADD_{}_*.root'.format(box))
    done_file_name = 'done_limits_{}_{}.txt'.format(model, box)
    with open(done_file_name, 'w') as done_file:
        for f in done_jobs:
            done_file.write(os.path.basename(f)+'\n')
    return done_file_name

def submit_box(model, box, no_sub=True):
    """
    Submit jobs for a single box.
    Do not resubmit done jobs.
    """
    script = 'python/RunMADDLimitJobs.py'
    out_dir = VERSION
    queue = '8nh'
    done_file_name = write_done_file(model, box)
    command = ['python', script, '--box', box, '--model', model,
            '--dir', out_dir, '--queue', queue, 
            '--done-file', done_file_name]
    if no_sub:
        command.append('--no-sub')
    do_command(command, False)

def submit(model, tag, sms, no_sub=True):
    """
    Submits the limit jobs for one SMS scan.
    model: string - name of model (without 'SMS-')
    tag: e.g. Razor2016_MoriondRereco
    sms: SMS object containing scan information
    no_sub: if True, do not submit jobs
    """
    for box in sms.boxes:
        print "Box {}".format(box)
        if sms.submodels is not None:
            for submodel in sms.submodels:
                print "Submitting jobs for {}".format(submodel)
                submit_box(submodel, box, no_sub)
        else:
            submit_box(model, box, no_sub)

def submit_combine(model, tag, sms, no_sub=True):
    """
    Submits jobs to combine limits for one SMS scan.
    Combines all relevant boxes for the chosen model.
    """
    script = 'python/CombineMADDLimits.py'
    out_dir = get_limit_dir(model)
    queue = '1nd'
    done_file_name = write_done_file(model, '_'.join(sms.boxes))
    command = ['python', script, '--model', model,
            '--dir', out_dir, '--queue', queue,
            '--done-file', done_file_name, '--boxes']
    command += sms.boxes
    if no_sub:
        command.append('--no-sub')
    do_command(command, False)

def aggregate(model, tag, sms, no_exec=True):
    """
    Makes a directory containing results from all submodels
    for a given SMS. Rename job files in the new directory
    to uniformize the output files.
    """
    if sms.submodels is None:
        print "No submodels to aggregate for {}".format(model)
        return
    out_dir = get_limit_dir(model)
    do_command(['mkdir', '-p', out_dir], no_exec)
    for submodel in sms.submodels:
        print "Dataset: {}".format(submodel)
        in_dir = get_limit_dir(submodel)
        in_files = glob.glob(in_dir
            +'/higgsCombineMADD_*.root')
        in_files += glob.glob(in_dir
            +'/RazorInclusiveMADD_*.*')
        for f in in_files:
            out_f = f.replace(submodel, model)
            if os.path.isfile(out_f):
                print "File {} already exists! Check for duplicate mass points."
            if f.endswith('.txt'):
                # change the name of the histogram file within the card
                sed_cmd = "sed 's/{}/{}/g' {} > {}".format(submodel, model,
                        f, out_f)
                print sed_cmd
                if not no_exec:
                    os.system(sed_cmd)
            else:
                do_command(['cp', f, out_f], no_exec)

def get_plot_dir():
    return 'PlotsSMS/plots/{}'.format(VERSION)

def get_config_dir():
    return 'PlotsSMS/config/{}'.format(VERSION)

def get(model, tag, sms, no_smooth=False, do_combined=False,
        no_exec=True):
    """
    Retrieves limit results.  Args are the same as submit()
    except:
    no_smooth: do not fill gaps between points
    do_combined: retrieve combined limit 
    """
    get_script = 'python/GetCombineMADD.py'
    contour_script = 'python/Get2DContour.py'
    out_dir = get_limit_dir(model)
    boxes = sms.boxes
    if do_combined:
        boxes = ['_'.join(sms.boxes)]
    for box in boxes:
        print "Box {}".format(box)
        command = ['python', get_script, '--box', box, '--model', model, 
                '--dir', out_dir]
        do_command(command, no_exec)
        
        command = ['python', contour_script, '--box', box, '--model', model, 
                '--dir', out_dir]
        if not sms.isGluino:
            command += ['--xsec-file', 'data/stop13TeV.txt']
        if no_smooth:
            command.append('--no-smooth')
        do_command(command, no_exec)

def get_box_str(box):
    if 'Lepton' not in box:
        return 'razor_0L'
    for piece in box.split('_'):
        if piece == 'DiJet' or piece == 'MultiJet':
            return 'razor_0L+1L'
    return 'razor_1L'

def write_plot_config(model, box, blind=True, preliminary=True, 
        lumi=35900, no_exec=True):
    """
    Writes a config file for making SMS plot.
    Returns the name of the config.
    """
    in_dir = get_limit_dir(model)
    name = '{}_{}'.format(model, box)
    results_file = '{}/{}_results.root'.format(
            in_dir, name)
    if blind:
        hist_type = 'Exp'
    else:
        hist_type = 'Obs'
    hist_name = 'xsecUL_{}_{}'.format(hist_type, name)
    text = "HISTOGRAM {} {}\n".format(results_file, hist_name)
    text += "EXPECTED {} Exp_{} ExpPlus_{} ExpMinus_{} kRed kOrange\n".format(
            results_file, name, name, name)
    text += "EXPECTED2 {} Exp_{} ExpPlus2_{} ExpMinus2_{} kGray+1 kGray\n".format(
            results_file, name, name, name)
    if not blind:
        text += "OBSERVED {} Obs_{} ObsPlus_{} ObsMinus_{} kBlack kBlue-9\n".format(
                results_file, name, name, name)
    text += "PRELIMINARY"
    if preliminary:
        text += " preliminary\n"
    else:
        text += " \n"
    text += "LUMI {}\n".format(lumi)
    text += "ENERGY 13\n"
    box_str = get_box_str(box)
    text += "BOXES {}\n".format(box_str)

    print "Plotting config:\n"
    print text
    out_dir = get_config_dir()
    out_f = '{}/{}.config'.format(out_dir, name)
    if blind:
        out_f = out_f.replace('.config', '_blinded.config')
    if preliminary:
        out_f = out_f.replace('.config', '_preliminary.config')
    if not no_exec:
        do_command(['mkdir', '-p', out_dir], no_exec)
        with open(out_f, 'w') as f:
            f.write(text)
        print "Plot config written to {}".format(out_f)
    return out_f

def plot(model, tag, sms, blind=True, preliminary=True, no_smooth=False,
        no_exec=True, do_combined=False):
    """
    Plots limit results.  Args are the same as submit() except for:
    blind: do not draw observed limit
    preliminary: write 'preliminary' on plots
    no_smooth: filename indicates that points were not smoothed.
    """
    plot_dir = get_plot_dir()
    do_command(['mkdir', '-p', plot_dir], no_exec)
    boxes = sms.boxes 
    if do_combined:
        boxes = ['_'.join(sms.boxes)]
    for box in boxes:
        print "Box: {}".format(box)
        config = write_plot_config(model, box, blind, preliminary, 
                no_exec=no_exec)
        name = '{}/{}_{}_'.format(plot_dir, model, box)
        if blind:
            name += 'Blinded_'
        if preliminary:
            name += 'Preliminary_'
        if no_smooth:
            name += 'NoSmooth_'
        script = 'PlotsSMS/python/makeSMSplots.py'
        command = ['python', script, config, name]
        do_command(command, no_exec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model')
    parser.add_argument('--tag', default='Razor2016_MoriondRereco')
    parser.add_argument('--no-exec', action='store_true')
    parser.add_argument('--no-sub', action='store_true',
            help='same as --no-exec')

    # actions
    parser.add_argument('--submit', action='store_true')
    parser.add_argument('--combined', action='store_true',
            help='do combination of boxes')
    parser.add_argument('--get', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--aggregate', action='store_true',
            help='combine datasets (for SMS split into multiple datasets')
    parser.add_argument('--finish', action='store_true',
            help='equivalent to --aggregate --get --plot')

    # configuration options
    parser.add_argument('--blind', action='store_true')
    parser.add_argument('--preliminary', action='store_true')
    parser.add_argument('--no-smooth', action='store_true',
            help='do not fill gaps between plotted points')

    args = parser.parse_args()
    no_exec = (args.no_exec or args.no_sub)

    # avoid looking at unblinded limits yet
    args.blind = True
    args.preliminary = True

    print "Model: {}".format(args.model)
    print "Analysis tag: {}".format(args.tag)
    try:
        sms = sms_models[args.model]
    except KeyError:
        sys.exit("Model {} is not implemented!".format(args.model))

    if args.submit:
        print "Submit limit jobs for {}\n".format(args.model)
        if args.combined:
            submit_combine(args.model, args.tag, sms, no_exec)
        else:
            submit(args.model, args.tag, sms, no_exec)
        sys.exit()

    if (args.aggregate or args.finish) and not args.combined:
        if sms.submodels is not None:
            print "Combine limit jobs from all datasets for {}".format(
                    args.model)
            aggregate(args.model, args.tag, sms, no_exec)
        else:
            print "Ignoring --aggregate, no submodels to aggregate"

    if args.get or args.finish:
        print "Get limit results for {}\n".format(args.model)
        get(args.model, args.tag, sms, no_smooth=args.no_smooth, 
                do_combined=args.combined, no_exec=no_exec)

    if args.plot or args.finish:
        print "Plot limit for {}\n".format(args.model)
        plot(args.model, args.tag, sms, blind=args.blind,
                preliminary=args.preliminary, no_smooth=args.no_smooth,
                do_combined=args.combined, no_exec=no_exec)
