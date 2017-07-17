"""
Module containing general utility methods.
"""
import ConfigParser, gzip, locale, os, random, subprocess, sys, time
from operator import itemgetter, attrgetter

def asList(value, delim=',') :
    """Generic method for creating a list from a given value.
    The value may be a string of delimiter-separated values, 
    a list or a set."""
    if type(value) == str :
        return value.split(delim)
    elif type(value) == list :
        return value
    elif type(value) == set :
        return list(value)
    else :
        raise ValueError('Expected a string, a list or a set; received %s' % type(value))


def commaFormat(d) :
    """Formats integer values using commas.  For example, 123456789 becomes '123,456,789'"""
    # Establish user's default their OS:
    locale.setlocale(locale.LC_ALL, '')
    return locale.format("%d", d, grouping=True)

def configMap(cfgFile):
    """Reads a configuration file and returns a dictionary of all sections, options and values."""
    result = {}
    config = ConfigParser.ConfigParser()
    # Activates case-sensitivity:
    config.optionxform = str

    try :
        config.read(cfgFile)
    except ConfigParser.ParsingError :
        raise ValueError('Invalid configuration file %s' % cfgFile)

    for sect in config.sections() :
        result[sect] = {}
        options = config.options(sect)
        for opt in options :
            try :
                result[sect][opt] = config.get(sect,opt)
            except :
                result[option] = None
    return result


def ezOpen(fileName) :
    """Allows clients to open files without regard for whether they're gzipped."""
    if not (os.path.exists(fileName) and os.path.isfile(fileName)):
        raise ValueError, 'file does not exist at %s' % fileName
    
    fileHandle  = gzip.GzipFile(fileName)
    gzippedFile = True
    try :
        line = fileHandle.readline()
        fileHandle.close()
        return gzip.GzipFile(fileName)
    except :
        return open(fileName)

def file2Dict(f,keyindex=0):
    result = {}
    fh=ezOpen(f)
    for line in fh:
        if line.startswith("#"):
            continue
        arr=line.rstrip().split("\t")
        result[arr[keyindex]]=line
    fh.close()
    return result


def filePrefix(f) :
    """Returns the filename prefix for a file.  For example:
       /my/dir/myfile.ext --> myfile"""
    head,tail     = os.path.split(f)
    prefix,suffix = os.path.splitext(tail)
    return prefix


def getAttribute(key, default, **args) :
    """Returns the value for the given key in the arguments dict, if found;
    otherwise returns default value."""
    return default if key not in args else args[key]

def getEnvironmentValue(name, default=None) :
    """Returns the value for the given environment variable name, if found;
    otherwise returns default value."""
    try :
        return os.environ[name]
    except KeyError :
        return default


def runCommand(s, **args) :
    """Announces a command runs it.."""
    logstream   = getAttribute('logstream', None, **args)
    debug       = getAttribute('debug', False, **args)
    exitOnError = getAttribute('exitOnError', True, **args)
    stderr      = getAttribute('stderr', None, **args)
    stdout      = getAttribute('stdout', None, **args)
    message     = '    ' + timeString('%s\n' % s)
    sys.stderr.write(message)
    if logstream : logstream.write(message)

    retcode = 0
    if not debug :
        if not stderr : streams.hideStderr()
        if not stdout : streams.hideStdout()
        retcode = subprocess.call(s, shell=True, stderr=stderr, stdout=stdout)
        if not stderr : streams.showStderr()
        if not stdout : streams.showStdout()

    if exitOnError and retcode < 0 :
        raise Exception('Error running command: returned %d signal\n%s' % (retcode, s))


def timeStamp(format_string='%Y%m%d%H%M%S') :
    """Returns a timestamp unique to the current second."""
    return time.strftime(format_string, time.localtime())

def timeString(s, format_string='%X', LF=False) :
    """Returns the input string with user-readable a timestamp prefix."""
    timestamp = time.strftime(format_string, time.localtime())
    result    =  '%s %s' % (timestamp, s)
    if LF : result += '\n'
    return result


def validateDir(path) :
    """Standard method for validating directory paths."""
    validateFile(path)
    if not os.path.isdir(path) :
        raise Exception("'%s' is not a directory; exiting." % path)

def validateFile(path) :
    """Standard method for validating file paths."""
    if not path :
        raise Exception("'%s' is not a valid file path; exiting." % path)

    if not os.path.exists(path) :
        raise Exception("File '%s' not found; exiting." % path)


def sortArr(arr,*args):
    sortarr=sorted(arr, key=itemgetter(*args))
    return sortarr

def revsortArr(arr,*args):
    sortarr=sorted(arr, key=itemgetter(*args),reverse=True)
    return sortarr


class ProgressIndicator(object) :
    """A simple progress indicator."""
    def __init__(self, increment, description='', verbose=True) :
        self.limit   = increment
        self.barlim  = int(self.limit/10)
        self.dotlim  = int(self.barlim/5)
        self.descr   = description
        self.started = False
        self.verbose = verbose
        self.ctr     = 0

    def count(self) :
        """Returns the current count."""
        return self.ctr

    def finish(self) :
        """Finishes the progress output by appending a newline,
        if anything has been written."""
        if self.started : sys.stderr.write('\n')
        self.started = False

    def reset(self) :
        """Resets the indicator to be used again."""
        self.ctr = 0
        self.finish()

    def update(self) :
        """Updates the indicator."""
        self.ctr += 1
        if not self.verbose : return
        if self.ctr % self.limit == 0 :
            sys.stderr.write('%s %s\n' % (commaFormat(self.ctr), self.descr))
        elif self.ctr % self.barlim == 0 :
            sys.stderr.write('|')
        elif self.ctr % self.dotlim == 0 :
            self.started = True
            sys.stderr.write('.')

