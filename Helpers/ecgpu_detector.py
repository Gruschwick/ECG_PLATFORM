import copy
import pdb
import numpy as np
from scipy import signal
from sklearn.preprocessing import normalize

from wfdb.processing.basic import get_filter_gain
from wfdb.processing.peaks import find_local_peaks
from wfdb.io.record import Record

class XQRS(object):
    """
    The QRS detector class for the XQRS algorithm. The `XQRS.Conf` 
    class is the configuration class that stores initial parameters 
    for the detection. The `XQRS.detect` method runs the detection algorithm.

    The process works as follows:

    - Load the signal and configuration parameters.
    - Bandpass filter the signal between 5 and 20 Hz, to get the
      filtered signal.
    - Apply moving wave integration (MWI) with a Ricker
      (Mexican hat) wavelet onto the filtered signal, and save the
      square of the integrated signal.
    - Conduct learning if specified, to initialize running
      parameters of noise and QRS amplitudes, the QRS detection
      threshold, and recent R-R intervals. If learning is unspecified
      or fails, use default parameters. See the docstring for the
      `_learn_init_params` method of this class for details.
    - Run the main detection. Iterate through the local maxima of
      the MWI signal. For each local maxima:

      - Check if it is a QRS complex. To be classified as a QRS,
        it must come after the refractory period, cross the QRS
        detection threshold, and not be classified as a T-wave
        if it comes close enough to the previous QRS. If
        successfully classified, update running detection
        threshold and heart rate parameters.
      - If not a QRS, classify it as a noise peak and update
        running parameters.
      - Before continuing to the next local maxima, if no QRS
        was detected within 1.66 times the recent R-R interval,
        perform backsearch QRS detection. This checks previous
        peaks using a lower QRS detection threshold.

    Attributes
    ----------
    sig : 1d ndarray
        The input ECG signal to apply the QRS detection on.
    fs : int, float
        The sampling frequency of the input signal.
    conf : XQRS.Conf object, optional
        The configuration object specifying signal configuration
        parameters. See the docstring of the XQRS.Conf class.

    Examples
    --------
    >>> import wfdb
    >>> from wfdb import processing

    >>> sig, fields = wfdb.rdsamp('sample-data/100', channels=[0])
    >>> xqrs = processing.XQRS(sig=sig[:,0], fs=fields['fs'])
    >>> xqrs.detect()

    >>> wfdb.plot_items(signal=sig, ann_samp=[xqrs.qrs_inds])

    """
    def __init__(self, sig, fs, conf=None):
        if sig.ndim != 1:
            raise ValueError('sig must be a 1d numpy array')
        self.sig = sig
        self.fs = fs
        self.sig_len = len(sig)
        self.conf = conf or XQRS.Conf()
        self._set_conf()


    class Conf(object):
        """
        Initial signal configuration object for this QRS detector.

        Attributes
        ----------
        hr_init : int, float, optional
            Initial heart rate in beats per minute. Used for calculating
            recent R-R intervals. 
        hr_max : int, float, optional
            Hard maximum heart rate between two beats, in beats per
            minute. Used for refractory period.
        hr_min : int, float, optional
            Hard minimum heart rate between two beats, in beats per
            minute. Used for calculating recent R-R intervals.
        qrs_width : int, float, optional
            Expected QRS width in seconds. Used for filter widths
            indirect refractory period.
        qrs_thr_init : int, float, optional
            Initial QRS detection threshold in mV. Use when learning
            is False, or learning fails.
        qrs_thr_min : int, float, string, optional
            Hard minimum detection threshold of QRS wave. Leave as 0
            for no minimum.
        ref_period : int, float, optional
            The QRS refractory period.
        t_inspect_period : int, float, optional
            The period below which a potential QRS complex is
            inspected to see if it is a T-wave.

        """
        def __init__(self, hr_init=75, hr_max=200, hr_min=25, qrs_width=0.1,
                     qrs_thr_init=0.13, qrs_thr_min=0, ref_period=0.2,
                     t_inspect_period=0.36):
            if hr_min < 0:
                raise ValueError("'hr_min' must be >= 0")

            if not hr_min < hr_init < hr_max:
                raise ValueError("'hr_min' < 'hr_init' < 'hr_max' must be True")

            if qrs_thr_init < qrs_thr_min:
                raise ValueError("qrs_thr_min must be <= qrs_thr_init")

            self.hr_init = hr_init
            self.hr_max = hr_max
            self.hr_min = hr_min
            self.qrs_width = qrs_width
            self.qrs_radius = self.qrs_width / 2
            self.qrs_thr_init = qrs_thr_init
            self.qrs_thr_min = qrs_thr_min
            self.ref_period = ref_period
            self.t_inspect_period = t_inspect_period


    def _set_conf(self):
        """
        Set configuration parameters from the Conf object into the detector
        object. Time values are converted to samples, and amplitude values 
        are in mV.

        Parameters
        ----------
        N/A

        Returns
        -------
        N/A

        """
        self.rr_init = 60 * self.fs / self.conf.hr_init
        self.rr_max = 60 * self.fs / self.conf.hr_min
        self.rr_min = 60 * self.fs / self.conf.hr_max

        # Note: if qrs_width is odd, qrs_width == qrs_radius*2 + 1
        self.qrs_width = int(self.conf.qrs_width * self.fs)
        self.qrs_radius = int(self.conf.qrs_radius * self.fs)

        self.qrs_thr_init = self.conf.qrs_thr_init
        self.qrs_thr_min = self.conf.qrs_thr_min

        self.ref_period = int(self.conf.ref_period * self.fs)
        self.t_inspect_period = int(self.conf.t_inspect_period * self.fs)


    def _bandpass(self, fc_low=5, fc_high=20):
        """
        Apply a bandpass filter onto the signal, and save the filtered
        signal.

        Parameters
        ----------
        fc_low : int, float
            The low frequency cutoff for the filter.
        fc_high : int, float
            The high frequency cutoff for the filter.

        Returns
        -------
        N/A

        """
        self.fc_low = fc_low
        self.fc_high = fc_high

        b, a = signal.butter(2, [float(fc_low) * 2 / self.fs,
                                 float(fc_high) * 2 / self.fs], 'pass')
        self.sig_f = signal.filtfilt(b, a, self.sig[self.sampfrom:self.sampto],
                                     axis=0)
        # Save the passband gain (x2 due to double filtering)
        self.filter_gain = get_filter_gain(b, a, np.mean([fc_low, fc_high]),
                                           self.fs) * 2


    def _mwi(self):
        """
        Apply moving wave integration (MWI) with a Ricker (Mexican hat)
        wavelet onto the filtered signal, and save the square of the
        integrated signal. The width of the hat is equal to the QRS width.
        After integration, find all local peaks in the MWI signal.

        Parameters
        ----------
        N/A

        Returns
        -------
        N/A

        """
        wavelet_filter = signal.ricker(self.qrs_width, 4)

        self.sig_i = signal.filtfilt(wavelet_filter, [1], self.sig_f,
                                     axis=0) ** 2

        # Save the MWI gain (x2 due to double filtering) and the total
        # gain from raw to MWI
        self.mwi_gain = get_filter_gain(wavelet_filter, [1],
                         np.mean([self.fc_low, self.fc_high]), self.fs) * 2
        self.transform_gain = self.filter_gain * self.mwi_gain
        self.peak_inds_i = find_local_peaks(self.sig_i, radius=self.qrs_radius)
        self.n_peaks_i = len(self.peak_inds_i)


    def _learn_init_params(self, n_calib_beats=8):
        """
        Find a number of consecutive beats and use them to initialize:
        - recent QRS amplitude
        - recent noise amplitude
        - recent R-R interval
        - QRS detection threshold

        The learning works as follows:
        - Find all local maxima (largest sample within `qrs_radius`
          samples) of the filtered signal.
        - Inspect the local maxima until `n_calib_beats` beats are
          found:
          - Calculate the cross-correlation between a Ricker wavelet of
            length `qrs_width`, and the filtered signal segment centered
            around the local maximum.
          - If the cross-correlation exceeds 0.6, classify it as a beat.
        - Use the beats to initialize the previously described
          parameters.
        - If the system fails to find enough beats, the default
          parameters will be used instead. See the docstring of
          `XQRS._set_default_init_params` for details.

        Parameters
        ----------
        n_calib_beats : int, optional
            Number of calibration beats to detect for learning

        Returns
        -------
        N/A

        """
        if self.verbose:
            print('Learning initial signal parameters...')

        last_qrs_ind = -self.rr_max
        qrs_inds = []
        qrs_amps = []
        noise_amps = []

        ricker_wavelet = signal.ricker(self.qrs_radius * 2, 4).reshape(-1,1)

        # Find the local peaks of the signal.
        peak_inds_f = find_local_peaks(self.sig_f, self.qrs_radius)

        # Peak numbers at least qrs_width away from signal boundaries
        peak_nums_r = np.where(peak_inds_f > self.qrs_width)[0]
        peak_nums_l = np.where(peak_inds_f <= self.sig_len - self.qrs_width)[0]

        # Skip if no peaks in range
        if (not peak_inds_f.size or not peak_nums_r.size
                                 or not peak_nums_l.size):
            if self.verbose:
                print('Failed to find %d beats during learning.'
                      % n_calib_beats)
            self._set_default_init_params()
            return

        # Go through the peaks and find QRS peaks and noise peaks.
        # only inspect peaks with at least qrs_radius around either side
        for peak_num in range(peak_nums_r[0], peak_nums_l[-1]):
            i = peak_inds_f[peak_num]
            # Calculate cross-correlation between the filtered signal
            # segment and a Ricker wavelet

            # Question: should the signal be squared? Case for inverse QRS
            # complexes
            sig_segment = normalize((self.sig_f[i - self.qrs_radius:
                                                i + self.qrs_radius]).reshape(-1, 1), axis=0)

            xcorr = np.correlate(sig_segment[:, 0], ricker_wavelet[:,0])

            # Classify as QRS if xcorr is large enough
            if xcorr > 0.6 and i-last_qrs_ind > self.rr_min:
                last_qrs_ind = i
                qrs_inds.append(i)
                qrs_amps.append(self.sig_i[i])
            else:
                noise_amps.append(self.sig_i[i])

            if len(qrs_inds) == n_calib_beats:
                break

        # Found enough calibration beats to initialize parameters
        if len(qrs_inds) == n_calib_beats:

            if self.verbose:
                print('Found %d beats during learning.' % n_calib_beats
                      + ' Initializing using learned parameters')

            # QRS amplitude is most important.
            qrs_amp = np.mean(qrs_amps)

            # Set noise amplitude if found
            if noise_amps:
                noise_amp = np.mean(noise_amps)
            else:
                # Set default of 1/10 of QRS amplitude
                noise_amp = qrs_amp / 10

            # Get R-R intervals of consecutive beats, if any.
            rr_intervals = np.diff(qrs_inds)
            rr_intervals = rr_intervals[rr_intervals < self.rr_max]
            if rr_intervals.any():
                rr_recent = np.mean(rr_intervals)
            else:
                rr_recent = self.rr_init

            # If an early QRS was detected, set last_qrs_ind so that it can be
            # picked up.
            last_qrs_ind = min(0, qrs_inds[0] - self.rr_min - 1)

            self._set_init_params(qrs_amp_recent=qrs_amp,
                                  noise_amp_recent=noise_amp,
                                  rr_recent=rr_recent,
                                  last_qrs_ind=last_qrs_ind)
            self.learned_init_params = True

        # Failed to find enough calibration beats. Use default values.
        else:
            if self.verbose:
                print('Failed to find %d beats during learning.'
                      % n_calib_beats)

            self._set_default_init_params()


    def _set_init_params(self, qrs_amp_recent, noise_amp_recent, rr_recent,
                         last_qrs_ind):
        """
        Set initial online parameters.

        Parameters
        ----------
        qrs_amp_recent : int, float
            The mean of the signal QRS amplitudes.
        noise_amp_recent : int, float
            The mean of the signal noise amplitudes.
        rr_recent : int
            The mean of the signal R-R interval values.
        last_qrs_ind : int
            The index of the signal's early QRS detected.

        Returns
        -------
        N/A

        """
        self.qrs_amp_recent = qrs_amp_recent
        self.noise_amp_recent = noise_amp_recent
        # What happens if qrs_thr is calculated to be less than the explicit
        # min threshold? Should print warning?
        self.qrs_thr = max(0.25*self.qrs_amp_recent
                           + 0.75*self.noise_amp_recent,
                           self.qrs_thr_min * self.transform_gain)
        self.rr_recent = rr_recent
        self.last_qrs_ind = last_qrs_ind

        # No QRS detected initially
        self.last_qrs_peak_num = None


    def _set_default_init_params(self):
        """
        Set initial running parameters using default values.

        The steady state equation is:
          `qrs_thr = 0.25*qrs_amp + 0.75*noise_amp`

        Estimate that QRS amp is 10x noise amp, giving:
          `qrs_thr = 0.325 * qrs_amp or 13/40 * qrs_amp`

        Parameters
        ----------
        N/A

        Returns
        -------
        N/A

        """
        if self.verbose:
            print('Initializing using default parameters')
        # Multiply the specified ECG thresholds by the filter and MWI gain
        # factors
        qrs_thr_init = self.qrs_thr_init * self.transform_gain
        qrs_thr_min = self.qrs_thr_min * self.transform_gain

        qrs_amp = 27/40 * qrs_thr_init
        noise_amp = qrs_amp / 10
        rr_recent = self.rr_init
        last_qrs_ind = 0

        self._set_init_params(qrs_amp_recent=qrs_amp,
                              noise_amp_recent=noise_amp,
                              rr_recent=rr_recent,
                              last_qrs_ind=last_qrs_ind)

        self.learned_init_params = False


    def _is_qrs(self, peak_num, backsearch=False):
        """
        Check whether a peak is a QRS complex. It is classified as QRS
        if it:
        - Comes after the refractory period.
        - Passes QRS threshold.
        - Is not a T-wave (check it if the peak is close to the previous QRS).

        Parameters
        ----------
        peak_num : int
            The peak number of the MWI signal to be inspected.
        backsearch: bool, optional
            Whether the peak is being inspected during backsearch.

        Returns
        -------
        bool
            Whether the peak is QRS (True) or not (False).

        """
        i = self.peak_inds_i[peak_num]
        if backsearch:
            qrs_thr = self.qrs_thr / 2
        else:
            qrs_thr = self.qrs_thr

        if (i-self.last_qrs_ind > self.ref_period
           and self.sig_i[i] > qrs_thr):
            if i-self.last_qrs_ind < self.t_inspect_period:
                if self._is_twave(peak_num):
                    return False
            return True

        return False


    def _update_qrs(self, peak_num, backsearch=False):
        """
        Update live QRS parameters. Adjust the recent R-R intervals and
        QRS amplitudes, and the QRS threshold.

        Parameters
        ----------
        peak_num : int
            The peak number of the MWI signal where the QRS is detected.
        backsearch: bool, optional
            Whether the QRS was found via backsearch.
        
        Returns
        -------
        N/A

        """
        i = self.peak_inds_i[peak_num]

        # Update recent R-R interval if the beat is consecutive (do this 
        # before updating self.last_qrs_ind)
        rr_new = i - self.last_qrs_ind
        if rr_new < self.rr_max:
            self.rr_recent = 0.875*self.rr_recent + 0.125*rr_new

        self.qrs_inds.append(i)
        self.last_qrs_ind = i
        # Peak number corresponding to last QRS
        self.last_qrs_peak_num = self.peak_num

        # QRS recent amplitude is adjusted twice as quickly if the peak
        # was found via backsearch
        if backsearch:
            self.backsearch_qrs_inds.append(i)
            self.qrs_amp_recent = (0.75*self.qrs_amp_recent
                                   + 0.25*self.sig_i[i])
        else:
            self.qrs_amp_recent = (0.875*self.qrs_amp_recent
                                   + 0.125*self.sig_i[i])

        self.qrs_thr = max((0.25*self.qrs_amp_recent
                            + 0.75*self.noise_amp_recent), self.qrs_thr_min)

        return


    def _is_twave(self, peak_num):
        """
        Check whether a segment is a T-wave. Compare the maximum gradient of
        the filtered signal segment with that of the previous QRS segment.

        Parameters
        ----------
        peak_num : int
            The peak number of the MWI signal where the QRS is detected.
        
        Returns
        -------
        bool
            Whether a segment is a T-wave (True) or not (False).

        """
        i = self.peak_inds_i[peak_num]

        # Due to initialization parameters, last_qrs_ind may be negative.
        # No way to check in this instance.
        if self.last_qrs_ind - self.qrs_radius < 0:
            return False

        # Get half the QRS width of the signal to the left.
        # Should this be squared?
        sig_segment = normalize((self.sig_f[i - self.qrs_radius:i]
                                 ).reshape(-1, 1), axis=0)
        last_qrs_segment = self.sig_f[self.last_qrs_ind - self.qrs_radius:
                                      self.last_qrs_ind]

        segment_slope = np.diff(sig_segment)
        last_qrs_slope = np.diff(last_qrs_segment)

        # Should we be using absolute values?
        if max(segment_slope) < 0.5*max(abs(last_qrs_slope)):
            return True
        else:
            return False


    def _update_noise(self, peak_num):
        """
        Update live noise parameters.

        Parameters
        ----------
        peak_num : int
            The peak number.

        Returns
        -------
        N/A

        """
        i = self.peak_inds_i[peak_num]
        self.noise_amp_recent = (0.875*self.noise_amp_recent
                                 + 0.125*self.sig_i[i])
        return

    def _require_backsearch(self):
        """
        Determine whether a backsearch should be performed on prior peaks.

        Parameters
        ----------
        N/A

        Returns
        -------
        bool
            Whether to require backsearch (True) or not (False).

        """
        if self.peak_num == self.n_peaks_i-1:
            # If we just return false, we may miss a chance to backsearch.
            # Update this?
            return False

        next_peak_ind = self.peak_inds_i[self.peak_num + 1]

        if next_peak_ind-self.last_qrs_ind > self.rr_recent*1.66:
            return True
        else:
            return False


    def _backsearch(self):
        """
        Inspect previous peaks from the last detected QRS peak (if any),
        using a lower threshold.

        Parameters
        ----------
        N/A

        Returns
        -------
        N/A

        """
        if self.last_qrs_peak_num is not None:
            for peak_num in range(self.last_qrs_peak_num + 1, self.peak_num + 1):
                if self._is_qrs(peak_num=peak_num, backsearch=True):
                    self._update_qrs(peak_num=peak_num, backsearch=True)
                # No need to update noise parameters if it was classified as
                # noise. It would have already been updated.

    def _run_detection(self):
        """
        Run the QRS detection after all signals and parameters have been
        configured and set.

        Parameters
        ----------
        N/A

        Returns
        -------
        N/A

        """
        if self.verbose:
            print('Running QRS detection...')

        # Detected QRS indices
        self.qrs_inds = []
        # QRS indices found via backsearch
        self.backsearch_qrs_inds = []

        # Iterate through MWI signal peak indices
        for self.peak_num in range(self.n_peaks_i):
            if self._is_qrs(self.peak_num):
                self._update_qrs(self.peak_num)
            else:
                self._update_noise(self.peak_num)

            # Before continuing to the next peak, do backsearch if
            # necessary
            if self._require_backsearch():
                self._backsearch()

        # Detected indices are relative to starting sample
        if self.qrs_inds:
            self.qrs_inds = np.array(self.qrs_inds) + self.sampfrom
        else:
            self.qrs_inds = np.array(self.qrs_inds)

        if self.verbose:
            print('QRS detection complete.')


    def detect(self, sampfrom=0, sampto='end', learn=True, verbose=True):
        """
        Detect QRS locations between two samples.

        Parameters
        ----------
        sampfrom : int, optional
            The starting sample number to run the detection on.
        sampto : int, optional
            The final sample number to run the detection on. Set as
            'end' to run on the entire signal.
        learn : bool, optional
            Whether to apply learning on the signal before running the
            main detection. If learning fails or is not conducted, the
            default configuration parameters will be used to initialize
            these variables. See the `XQRS._learn_init_params` docstring
            for details.
        verbose : bool, optional
            Whether to display the stages and outcomes of the detection
            process.

        Returns
        -------
        N/A

        """
        if sampfrom < 0:
            raise ValueError("'sampfrom' cannot be negative")
        self.sampfrom = sampfrom

        if sampto == 'end':
            sampto = self.sig_len
        elif sampto > self.sig_len:
            raise ValueError("'sampto' cannot exceed the signal length")
        self.sampto = sampto
        self.verbose = verbose

        # Don't attempt to run on a flat signal
        if np.max(self.sig) == np.min(self.sig):
            self.qrs_inds = np.empty(0)
            if self.verbose:
                print('Flat signal. Detection skipped.')
            return

        # Get/set signal configuration fields from Conf object
        self._set_conf()
        # Bandpass filter the signal
        self._bandpass()
        # Compute moving wave integration of filtered signal
        self._mwi()

        # Initialize the running parameters
        if learn:
            self._learn_init_params()
        else:
            self._set_default_init_params()

        # Run the detection
        self._run_detection()


def xqrs_detect(sig, fs, sampfrom=0, sampto='end', conf=None,
                learn=True, verbose=True):
    """
    Run the 'xqrs' QRS detection algorithm on a signal. See the
    docstring of the XQRS class for algorithm details.

    Parameters
    ----------
    sig : ndarray
        The input ECG signal to apply the QRS detection on.
    fs : int, float
        The sampling frequency of the input signal.
    sampfrom : int, optional
        The starting sample number to run the detection on.
    sampto : str
        The final sample number to run the detection on. Set as 'end' to
        run on the entire signal.
    conf : XQRS.Conf object, optional
        The configuration object specifying signal configuration
        parameters. See the docstring of the XQRS.Conf class.
    learn : bool, optional
        Whether to apply learning on the signal before running the main
        detection. If learning fails or is not conducted, the default
        configuration parameters will be used to initialize these
        variables.
    verbose : bool, optional
        Whether to display the stages and outcomes of the detection
        process.

    Returns
    -------
    qrs_inds : ndarray
        The indices of the detected QRS complexes.

    Examples
    --------
    >>> import wfdb
    >>> from wfdb import processing

    >>> sig, fields = wfdb.rdsamp('sample-data/100', channels=[0])
    >>> qrs_inds = processing.xqrs_detect(sig=sig[:,0], fs=fields['fs'])

    """
    xqrs = XQRS(sig=sig, fs=fs, conf=conf)
    xqrs.detect(sampfrom=sampfrom, sampto=sampto, verbose=verbose)
    return xqrs.qrs_inds


def time_to_sample_number(seconds, frequency):
    """
    Convert time to sample number.

    Parameters
    ----------
    seconds : int, float
        The input time in seconds.
    frequency : int, float
        The input frequency.

    Returns
    -------
    float
        The converted sample number.

    """
    return seconds * frequency + 0.5
