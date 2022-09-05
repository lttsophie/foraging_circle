from __future__ import absolute_import, division

import os  # handy system and path functions
import random
from math import sin, cos, radians, sqrt

from psychopy import gui, visual, core, data, event, logging

# Ensure that relative paths start from the same directory as this script

_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2021.2.3'
expName = 'foraging_circle'  # from the Builder filename that created this script
expInfo = {'participant': '', 'age': '', 'gender': ''}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if not dlg.OK:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psy_exp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='1.1',
                                 extraInfo=expInfo, runtimeInfo=None,
                                 savePickle=True, saveWideText=True,
                                 dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename + '.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

# Window size
win_height = 1920
win_width = 1080

# Set up the Window
win = visual.Window(
    size=(win_height, win_width), fullscr=False, screen=0,
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-1, -1, -1], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    units='pix')

# Initialize components for Routine "Instruction"
InstructionClock = core.Clock()
mouse = event.Mouse(visible=True, win=win)

# Create text objects
instr_text_en = 'Please read the instructions for this experiment that the experimenter provided to you. ' \
                'If you have any questions, ask the experimenter. If everything is clear you can start ' \
                'the training session by pressing SPACE button. ' \
                'Try to search for objects as quickly and accurately as possible.'

end_train_text_en = 'The training session is over, if you have any questions, ask the experimenter. ' \
                    'To start the main session, press the SPACE button. ' \
                    'Try to search for objects as quickly and accurately as possible.'

instruction = visual.TextStim(win=win, name='text',
                              text=instr_text_en,
                              font='Open Sans',
                              pos=(0, 0), height=30, wrapWidth=1500, ori=0.0,
                              color='white', colorSpace='rgb', opacity=None,
                              languageStyle='LTR',
                              depth=0.0)

end_of_train = visual.TextStim(win=win, name='text',
                               text=end_train_text_en,
                               font='Open Sans',
                               pos=(0, 0), height=40, wrapWidth=1500, ori=0.0,
                               color='white', colorSpace='rgb', opacity=None,
                               languageStyle='LTR',
                               depth=0.0)

# Create default stimuli
green_cc = visual.Circle(win,
                         radius=20,
                         edges=32,
                         pos=(0, 0),
                         units='pix',
                         lineWidth=1.5,
                         fillColor='Green',
                         colorSpace='rgb')
# make several dictionaries for each object to know what is the name of the object
green_dict_cc = {'name': 'green_cc', 'figure': green_cc}

red_cc = visual.Circle(win,
                       radius=20,
                       edges=32,
                       pos=(0, 0),
                       units='pix',
                       lineWidth=1.5,
                       fillColor='Red',
                       colorSpace='rgb')
red_dict_cc = {'name': 'red_cc', 'figure': red_cc}

yellow_cc = visual.Circle(win,
                          radius=20,
                          edges=32,
                          pos=(0, 0),
                          units='pix',
                          lineWidth=1.5,
                          fillColor='Yellow',
                          colorSpace='rgb')
yellow_dict_cc = {'name': 'yellow_cc', 'figure': yellow_cc}

blue_cc = visual.Circle(win,
                        radius=20,
                        edges=32,
                        pos=(0, 0),
                        units='pix',
                        lineWidth=1.5,
                        fillColor='Blue',
                        colorSpace='rgb')
blue_dict_cc = {'name': 'blue_cc', 'figure': blue_cc}

green_sq = visual.Rect(win,
                       width=40,
                       height=40,
                       pos=(0, 0),
                       units='pix',
                       lineWidth=1.5,
                       fillColor='Green',
                       colorSpace='rgb')
green_dict_sq = {'name': 'green_sq', 'figure': green_sq}

red_sq = visual.Rect(win,
                     width=40,
                     height=40,
                     pos=(0, 0),
                     units='pix',
                     lineWidth=1.5,
                     fillColor='Red',
                     colorSpace='rgb')
red_dict_sq = {'name': 'red_sq', 'figure': red_sq}

# make list of dictionaries containing object names and objects itself [{name: red, figure: obj}]
stim_set_feat = [green_dict_cc, red_dict_cc, yellow_dict_cc, blue_dict_cc]
stim_set_conj = [green_dict_cc, red_dict_cc, green_dict_sq, red_dict_sq]

available_radius_for_click = 30


def make_lst_stims_in_circle(targ_list, condition):
    stim_set = [targ_list[0], targ_list[1]]
    rand_stims = random.choices(condition, k=4)
    stim_set = stim_set + rand_stims
    # return list of dictionaries (6 in total for each stimulus)
    return stim_set


def make_circle_coord(x, y, R):
    coord_name = []
    coord_value = []
    coord_value_jitter = []
    r = 20
    edge_x_1 = 960 - r
    edge_x_2 = -960 + r
    edge_y_1 = 540 - r
    edge_y_2 = -540 + r

    # Positioning in circle (e.g 'left', 'right', 'upper left' always same for all patches)
    # Positioning: {stim_0:right}, {stim1:upper right}, {stim2:upper left},
    # {stim3:left}, {stim4:lower left}, {stim5: lower right}
    for coord in range(0, 6):
        x_n = R * cos(radians(60 * coord)) + x
        y_n = R * sin(radians(60 * coord)) + y
        if edge_x_1 < x_n or edge_x_2 > x_n or edge_y_1 < y_n or y_n < edge_y_2:
            return False, False
        x_j = x_n + random.uniform(-5.0, 5.0)
        y_j = y_n + random.uniform(-5.0, 5.0)
        coord_name.append('x' + str(coord))
        coord_value.append(x_n)
        coord_value_jitter.append(x_j)
        coord_name.append('y' + str(coord))
        coord_value.append(y_n)
        coord_value_jitter.append(y_j)
    coords = zip(coord_name, coord_value)
    jit_coords = zip(coord_name, coord_value_jitter)
    circle_coord = dict(coords)
    circle_coord_jit = dict(jit_coords)
    # return dictionary with x, y values from 0 to 5 (e.g. {x0: num, y0: num, x1: num ...})
    return circle_coord, circle_coord_jit


def feature_conjunction_check(targ):
    if targ == [red_dict_cc, green_dict_cc]:
        thisExp.addData('condition', 'feat')
    elif targ == [yellow_dict_cc, blue_dict_cc]:
        thisExp.addData('condition', 'feat')
    if targ == [red_dict_cc, green_dict_sq]:
        thisExp.addData('condition', 'conj')
    if targ == [red_dict_sq, green_dict_cc]:
        thisExp.addData('condition', 'conj')


def quit_exp():
    keys = event.getKeys(keyList=['escape'])  # collect key presses after first flip
    if keys:
        if 'escape' in keys:
            thisExp.nextEntry()  # if someone wants to escape the experiment
            thisExp.addData('trial_n', 'BLEEP BLOOP SOMEONE PRESSED ESC')
            thisExp.saveAsWideText(filename + '.csv', delim='auto')
            core.quit()


def how_many_targets(stim_set, targets):
    count = 0
    for stim in range(0, 5):
        if stim_set[stim] == targets[0] or stim_set[stim] == targets[1]:
            count += 1
    return count


def make_stim_set(feat_t_rg, s_feat, feat_t_yb, s_conj, conj_t_rcc_gsq, conj_t_rsq_gcc):
    if curr_block == feat_t_rg:
        stim_set = make_lst_stims_in_circle(feat_t_rg, s_feat)
        random.shuffle(stim_set)
        return stim_set, curr_block
    elif curr_block == feat_t_yb:
        stim_set = make_lst_stims_in_circle(feat_t_yb, s_feat)
        random.shuffle(stim_set)
        return stim_set, curr_block
    elif curr_block == conj_t_rcc_gsq:
        stim_set = make_lst_stims_in_circle(conj_t_rcc_gsq, s_conj)
        random.shuffle(stim_set)
        return stim_set, curr_block
    elif curr_block == conj_t_rsq_gcc:
        stim_set = make_lst_stims_in_circle(conj_t_rsq_gcc, s_conj)
        random.shuffle(stim_set)
        return stim_set, curr_block


trial_break_text = visual.TextStim(win=win, name='text',
                                   text='Wrong answer',
                                   font='Open Sans',
                                   pos=(0, 0), height=60, wrapWidth=None, ori=0.0,
                                   color='white', colorSpace='rgb', opacity=None,
                                   languageStyle='LTR',
                                   depth=0.0)


def trial_procedure(stim_set, circle_coord_jit, circle_coord, targets, edge):
    global central_x, central_y, selected_stim, not_correct_answer
    # Record of all already known data
    thisExp.addData('block_n', block_n)
    thisExp.addData('trial_n', trial_n)
    thisExp.addData('training_trials', training_trials)
    thisExp.addData('target1', targets[0].get('name'))
    thisExp.addData('target2', targets[1].get('name'))
    if edge:
        thisExp.addData('back_to_center', 'TRUE')
    else:
        thisExp.addData('back_to_center', 'FALSE')
    # Count and record how many targets in patch
    thisExp.addData('num_of_targets', how_many_targets(stim_set, targets))
    feature_conjunction_check(targets)
    # creating and recording all stimuli and their positions
    for stim_n in range(0, 6):
        stim_pos = [circle_coord_jit.get(('x' + str(stim_n))), circle_coord_jit.get(('y' + str(stim_n)))]
        stim_set[stim_n].get('figure').pos = stim_pos
        thisExp.addData('stim' + str(stim_n), stim_set[stim_n].get('name'))
        thisExp.addData('x_' + str(stim_n), stim_pos[0])
        thisExp.addData('y_' + str(stim_n), stim_pos[1])

    # wait for click
    while mouse.getPressed()[0] == 0:
        for stim in range(0, 6):
            stim_pos = [circle_coord_jit.get(('x' + str(stim))), circle_coord_jit.get(('y' + str(stim)))]
            stim_set[stim].get('figure').pos = stim_pos
            stim_set[stim].get('figure').draw()
        quit_exp()
        win.flip(clearBuffer=True)

    # collecting mouse response and record Reaction time
    mouse_resp_pos = mouse.getPos()
    resp_x = mouse_resp_pos[0]
    resp_y = mouse_resp_pos[1]
    click_time = timer.getTime()
    thisExp.addData('x_resp', resp_x)
    thisExp.addData('y_resp', resp_y)
    thisExp.addData('click_time', click_time)

    # infinitely asking whether there was a click
    while mouse.getPressed()[0] != 0:
        for stim in range(0, 6):
            stim_pos = [circle_coord_jit.get(('x' + str(stim))), circle_coord_jit.get(('y' + str(stim)))]
            stim_set[stim].get('figure').pos = stim_pos
            stim_set[stim].get('figure').draw()
        win.flip(clearBuffer=True)

    # set counters to check if click was on object and if target was selected
    is_click_valid = False
    correct_answer = False

    # check if click is valid by calculating radius from center of the object to mouse response position
    for pos in range(0, 6):
        x = circle_coord.get(('x' + str(pos)))
        y = circle_coord.get(('y' + str(pos)))
        x_j = circle_coord_jit.get(('x' + str(pos)))
        y_j = circle_coord_jit.get(('y' + str(pos)))
        r = sqrt((x_j - resp_x) ** 2 + (y_j - resp_y) ** 2)  # Find out how far the click from center of object

        # compare distance from click and available distance
        # (considering the actual radius of the object + some additional distance)
        if r < available_radius_for_click:
            # set new phantom center for new coordinates creation
            central_x = x
            central_y = y
            selected_stim_name = stim_set[pos].get('name')
            # check if target was selected
            if selected_stim_name == targets[0].get('name') or selected_stim_name == targets[1].get('name'):
                correct_answer = True
                selected_time = timer.getTime()
                thisExp.addData('direction', str(pos))
                # check if there was a switch
                if selected_stim_name == selected_stim:
                    thisExp.addData('switch', 0)
                else:
                    thisExp.addData('switch', 1)
                selected_stim = selected_stim_name
                thisExp.addData('selected_stim', selected_stim_name)
                thisExp.addData('x_targ', x_j)
                thisExp.addData('y_targ', y_j)
                thisExp.addData('targ_selected', 'TRUE')

            else:
                not_correct_answer += 1
                thisExp.addData('targ_selected', 'FALSE')
                thisExp.addData('selected_stim', selected_stim_name)
                thisExp.addData('x_targ', x_j)
                thisExp.addData('y_targ', y_j)

            is_click_valid = True

    if is_click_valid:
        thisExp.addData('valid_click', 'TRUE')
    else:
        thisExp.addData('valid_click', 'FALSE')

    thisExp.nextEntry()  # go to the next raw in data file

    # if click was not in any object range run the function again
    if not is_click_valid:
        correct_answer = trial_procedure(stim_set, circle_coord_jit, circle_coord, targets, edge)
    return correct_answer


def trials(t_1, t_2, targ_rg, targ_yb, targ_rcc_gsq, targ_rsq_gcc, conj_stim, feat_stim, n_trials):
    global central_x, central_y, trial_n
    for trial in range(0, n_trials):
        # -Start of the trial screen with stimuli- #
        trial_n = trial
        central_x = 0
        central_y = 0
        t_1.pos = (-40, 0)
        t_2.pos = (40, 0)
        t_1.draw()
        t_2.draw()
        win.flip(clearBuffer=True)
        core.wait(0.5)
        timer.reset()

        # loop to show a patch of 6 objects (one screen with 6 stimuli)
        for stim_screen in range(0, 15):
            timer.reset()  # Beginning of the time recording
            stim_set, block_targets = make_stim_set(targ_rg, feat_stim,
                                                    targ_yb, conj_stim,
                                                    targ_rcc_gsq,
                                                    targ_rsq_gcc)

            circle_coord, circle_coord_jitter = make_circle_coord(central_x, central_y, 90)
            edge_reached = False
            # check if the edge was reached
            if not circle_coord:
                edge_reached = True
                central_x = 0
                central_y = 0
                circle_coord, circle_coord_jitter = make_circle_coord(central_x, central_y, 90)
                # if participant reached the edge and patch of objects moved to center

            check = trial_procedure(stim_set, circle_coord_jitter, circle_coord, block_targets, edge_reached)

            # if not target selected
            if not check:
                break_text = visual.TextStim(win=win, name='text',
                                             text='Wrong answer',
                                             font='Open Sans',
                                             pos=(0, 0), height=60, wrapWidth=None, ori=0.0,
                                             color='white', colorSpace='rgb', opacity=None,
                                             languageStyle='LTR',
                                             depth=0.0)
                break_text.draw()
                win.flip()
                core.wait(1.0)
                break


# --Make 4 blocks for targets-- #
feat_targ_rg = [red_dict_cc, green_dict_cc]
feat_targ_yb = [yellow_dict_cc, blue_dict_cc]
conj_targ_rcc_gsq = [red_dict_cc, green_dict_sq]
conj_targ_rsq_gcc = [red_dict_sq, green_dict_cc]

# Make a list of all blocks
stim_feat = [red_dict_cc, green_dict_cc, yellow_dict_cc, blue_dict_cc]
stim_conj = [red_dict_cc, green_dict_sq, red_dict_sq, green_dict_cc]
all_blocks = [feat_targ_rg, feat_targ_yb, conj_targ_rcc_gsq, conj_targ_rsq_gcc]

# ----START OF THE EXPERIMENT---- #

# --Show instruction--#
instruction.draw()
win.flip()
event.waitKeys(keyList=['space'])

# Create global counters
central_x = 0
central_y = 0
not_correct_answer = 0
# Create timer to record each selection time
timer = core.Clock()
# Create timer to record overall time for the experiment
time_for_participant = core.Clock()
block_n = 0
trial_n = 0
selected_stim = ''

# --Start training trials-- #
training_trials = 1
for block in range(0, 4):
    # Create text message for the beginning of each block
    block_num_text = visual.TextStim(win=win, name='text',
                                     text='Block #' + str(block + 1),
                                     font='Open Sans',
                                     pos=(0, 150), height=50, wrapWidth=None, ori=0.0,
                                     color='white', colorSpace='rgb', opacity=None,
                                     languageStyle='LTR',
                                     depth=0.0)
    block_key_text = visual.TextStim(win=win, name='text',
                                     text='Press SPACE to continue',
                                     font='Open Sans',
                                     pos=(0, -200), height=50, wrapWidth=None, ori=0.0,
                                     color='white', colorSpace='rgb', opacity=None,
                                     languageStyle='LTR',
                                     depth=0.0)
    # -Show targets for this block- #
    curr_block = random.choice(all_blocks)
    all_blocks.remove(curr_block)
    block_n = block
    # show targets for this block
    target1 = curr_block[0].get('figure')
    target2 = curr_block[1].get('figure')
    target1.pos = (-40, 0)
    target2.pos = (40, 0)
    target1.draw()
    target2.draw()
    block_num_text.draw()
    block_key_text.draw()
    win.flip(clearBuffer=True)
    event.waitKeys(keyList='space')  # waiting for a key press in the beginning of each block

    trials(target1, target2, feat_targ_rg,
           feat_targ_yb, conj_targ_rcc_gsq,
           conj_targ_rsq_gcc, stim_conj, stim_feat, 3)

training_trials = 0

# End of the training screen
end_of_train.draw()
win.flip()
event.waitKeys(keyList=['space'])

# Create new blocks
stim_feat = [red_dict_cc, green_dict_cc, yellow_dict_cc, blue_dict_cc]
stim_conj = [red_dict_cc, green_dict_sq, red_dict_sq, green_dict_cc]
all_blocks = [feat_targ_rg, feat_targ_yb, conj_targ_rcc_gsq, conj_targ_rsq_gcc]

time_for_participant.reset()
not_correct_answer = 0
# ---Start of the Experimental trials--- #
for block in range(0, 4):
    # Create text message for the beginning of each block
    block_num_text = visual.TextStim(win=win, name='text',
                                     text='Block #' + str(block + 1),
                                     font='Open Sans',
                                     pos=(0, 150), height=50, wrapWidth=None, ori=0.0,
                                     color='white', colorSpace='rgb', opacity=None,
                                     languageStyle='LTR',
                                     depth=0.0)
    block_key_text = visual.TextStim(win=win, name='text',
                                     text='Press SPACE to continue',
                                     font='Open Sans',
                                     pos=(0, -200), height=50, wrapWidth=None, ori=0.0,
                                     color='white', colorSpace='rgb', opacity=None,
                                     languageStyle='LTR',
                                     depth=0.0)
    # -Show targets for this block- #
    curr_block = random.choice(all_blocks)
    all_blocks.remove(curr_block)
    block_n = block
    target1 = curr_block[0].get('figure')
    target2 = curr_block[1].get('figure')
    target1.pos = (-40, 0)
    target2.pos = (40, 0)
    block_num_text.draw()
    block_key_text.draw()
    target1.draw()
    target2.draw()
    win.flip(clearBuffer=True)
    event.waitKeys(keyList='space')
    trials(target1, target2, feat_targ_rg,
           feat_targ_yb, conj_targ_rcc_gsq,
           conj_targ_rsq_gcc, stim_conj, stim_feat, 32)

# save the output data
thisExp.saveAsWideText(filename + '.csv', delim='auto')

# Show thw additional info for participant
num_of_trials = 32 * 4
result_time_for_participant = time_for_participant.getTime()
correct_answer_percentage = 100 - (not_correct_answer / num_of_trials * 100)
end_of_exp_result_time_text = visual.TextStim(win=win, name='text',
                                              text='You completed the experiment in ' + str(
                                                  "{:.2f}".format(result_time_for_participant / 60)) + ' min',
                                              font='Open Sans',
                                              pos=(0, 20), height=30, wrapWidth=None, ori=0.0,
                                              color='white', colorSpace='rgb', opacity=None,
                                              languageStyle='LTR',
                                              depth=0.0)
end_of_exp_correct_answer_percentage_text = visual.TextStim(win=win, name='text',
                                                            text='Total percentage of correctly finished trials: ' +
                                                                 str(correct_answer_percentage) + '%',
                                                            font='Open Sans',
                                                            pos=(0, -60), height=30, wrapWidth=None, ori=0.0,
                                                            color='white', colorSpace='rgb', opacity=None,
                                                            languageStyle='LTR',
                                                            depth=0.0)
end_of_exp_text = visual.TextStim(win=win, name='text',
                                  text='Thank you for your participation! Press ESCAPE to finish.',
                                  font='Open Sans',
                                  pos=(0, -140), height=30, wrapWidth=None, ori=0.0,
                                  color='white', colorSpace='rgb', opacity=None,
                                  languageStyle='LTR',
                                  depth=0.0)
end_of_exp_result_time_text.draw()
end_of_exp_correct_answer_percentage_text.draw()
end_of_exp_text.draw()
win.flip()
end_key = event.waitKeys(keyList='escape')
