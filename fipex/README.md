Parsing des scripts FIPEX du VKI

    % make
    cc -Wall -Wextra -o su.o -c su.c
    cc -o su su.o

    % ./su < vki-script.bin 
        script length: 67 bytes
     start time (UTC): Wed Jan  1 12:00:00 2014
          repeat time: 3600s
        command count: 10
    CMD_ID: 0x0f, LEN: 0, XOR: 0x0f (match), DELAY: 60s
    CMD_ID: 0x0b, LEN: 0, XOR: 0x0b (match), DELAY: 60s
    CMD_ID: 0x11, LEN: 3, DATA: [ 0x04, 0x01, 0x00 ], XOR: 0x17 (match), DELAY: none
    CMD_ID: 0x11, LEN: 3, DATA: [ 0x05, 0x10, 0x0a ], XOR: 0x0d (match), DELAY: none
    CMD_ID: 0x11, LEN: 3, DATA: [ 0x02, 0xc8, 0x00 ], XOR: 0xd8 (match), DELAY: none
    CMD_ID: 0x0c, LEN: 0, XOR: 0x0c (match), DELAY: 300s
    CMD_ID: 0x20, LEN: 0, XOR: 0x20 (match), DELAY: none
    CMD_ID: 0x21, LEN: 0, XOR: 0x21 (match), DELAY: none
    CMD_ID: 0xf0, LEN: 0, XOR: 0xf0 (match), DELAY: none
    CMD_ID: 0xff
