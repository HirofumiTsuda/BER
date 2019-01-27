These files can be used to evaliate Bit Error Rate (BER) in Asynchronous CDMA systems.
This repository consists of the following files:
-- main functions --
・main.cpp
・db.cpp
-- library files --
・code.cpp
・code.hpp
・correlation.cpp
・correlation.hpp

1. main.cpp
With this file, we can invertigate BER in the following situation:
・SNR (the power of channel Gaussian Noise) is fixed. 
・The number of the users varies.

-- parameters (defined in #define)--
CHIP : the number of components in sequences. This parameter is commonly used in all the users.
BIT : the number of bits sent to the receiver.
MIN_USER_NUM : the minimum number of the users in the simulation.
MAX_USER_NUM : the maximum number of the users in the simulation.
INTERVAL : In this simulation, integral is calculated with trapezoidal approximation. This parameter denotes the width of each trapezoid. The calculation gets more precisely as this parameter gets larger.
LATE : the maximum delay-time (delay time is chosen randomly).
CHANNEL_LENGTH : the length of channel (This parameter must not be changed)
COV_COUNT : not used
USER_PLUS : the added number of users in each iteration.
MAX_FADING_NUM : not used (please do not change)
MAX_DURATION :  not used (please do not change)
FADING_POWER :  not used (please do not change)
SNR : Signal to Noise Ratio [db]
TRIAL_NUMBER : the number of iterations

-- how to use --
First, please set the parameters described above. Then, please choose the sequence.
This file requires the input file ``pnum.dat''.
The output file is defined in row-number 80. Please set your desired directory. 
Below the setBit (row number 123), there are some kinds of sequences. Please choose one.
These sequences are defined in ``code.cpp''.

-- Compile --
In this file, OpenMP technique is used. Further, correlation.o and code.o are required.
If you use Makefile, please type ``make''.

2. db.cpp
With this file, we can invertigate BER in the following situation:
・SNR (the power of channel Gaussian Noise) varies. 
・The number of the users is fixed.
Please notice that the above situation is different from one considered in main.cpp.

-- parameters (defined in #define)--
CHIP : the number of components in sequences. This parameter is commonly used in all the users.
BIT : the number of bits sent to the receiver.
USER_NUM : the number of the users in the simulation.
INTERVAL : In this simulation, integral is calculated with trapezoidal approximation. This parameter denotes the width of each trapezoid. The calculation gets more precisely as this parameter gets larger.
INIT_DB : initial SNR [db]
MAX_DB : the maximum SNR [db]
DB_PLUS : the added SNR [dB] in each iteration.
LATE : the maximum delay-time (delay time is chosen randomly).
CHANNEL_LENGTH : the length of channel (This parameter must not be changed)
COV_COUNT : not used
MAX_FADING_NUM : the maximum number of fading (delayed) signals.
MAX_DURATION :  not used (please do not change)
FADING_POWER :  the power of the fading (delayed) signals.
TRIAL_NUMBER : the number of iterations

-- how to use --
First, please set the parameters described above. Then, please choose the sequence.
This file requires the input file ``pnum.dat''.
The output file is defined in row-number 81. Please set your desired directory. 
Below the setBit (row number 128), there are some kinds of sequences. Please choose one.
These sequences are defined in ``code.cpp''.

-- Compile --
In this file, OpenMP technique is used. Further, correlation.o and code.o are required.
If you use Makefile, please type ``make''.

3. code.cpp
In this file, the sequences for asynchronous CDMA systems are defined.
If you would like to add sequences, please define and add your sequence to this file. Further, please add your function to ``code.hpp''.

4. correlation.cpp
Do not change.