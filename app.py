# ! /usr/bin/env python3
import pybedtools

from remus.remus import app
#pybedtools.debug_mode(True)

if __name__ == '__main__':
    app.run(host='0.0.0.0',
            debug=True)
