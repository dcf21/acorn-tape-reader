# acorn_tape_parse.py

This Python script extracts binary files from WAV recordings of audio cassette tapes recorded by 8-bit Acorn computers, including the BBC Micro, Acorn Electron and BBC Master computers. It can also produce UEF tape images for use in emulators such as BeebEm or JSBeeb.

This script was used by the author to recover all the Acorn tapes archived on the website <https://files.dcford.org.uk/>.

By default, this script simply converts a WAV recording into UEF format, and exports all the files to a specified output directory, together with a textual summary of the metadata associated with each file (load address, etc.). If a more sophisticated export is required, it is simple to call the <WavAcornFileSearch> class from an external script to perform other actions on the files found.

## Usage

* Any bit rate is supported, but >= 44.1kHz is recommended. Both mono and stereo recordings are accepted, but stereo is recommended, and the best channel will automatically be selected -- very often one channel is (much) less noisy than the other.

* This script currently assumes 1200 baud, as used by almost all software. It would probably be simple to add support for the BBC Micro's 300 baud setting, but the author has no suitable test data.

* This script does not currently support Acorn Atom tapes (which use a different header format). It would probably be simple to adapt it to do so, but the author has no suitable test data.

* This script is quite sensitive to low-frequency noise. You may be able to recover more files if you use audio-editing software (e.g. Adobe Audition or Audacity) to pass the input audio through a ~100-Hz high-pass filter before calling this script.

## Command-line syntax

./acorn_tape_parse.py [-h] [--input INPUT_FILENAME] [--output OUTPUT_DIR] [--uef FILENAME] [--debug] [--relative_speed RELATIVE_SPEED]

### Options:

| Switch                          | Meaning                                                                               |
|---------------------------------|---------------------------------------------------------------------------------------|
| -h, --help                      | Show help message and exit                                                            |
| --input INPUT_FILENAME          | Input WAV file to process                                                             |
| --output OUTPUT_DIR             | Directory in which to put the extracted files                                         |
| --uef FILENAME                  | Filename for output UEF tape image (for use in emulators such as BeebEm / JSBeeb)     |
| --debug                         | Show full debugging output                                                            |
| --relative_speed RELATIVE_SPEED | Assume tape runs at given fraction of normal speed (0.95 speeds up tape by 5 percent) |


## License

This code is distributed under the Gnu General Public License V3. It is (C) Dominic Ford 2022.

## Author

Dominic Ford - <https://dcford.org.uk/>
