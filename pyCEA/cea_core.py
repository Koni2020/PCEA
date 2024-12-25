"""
This is core utils in pyCEA package.
Created on 08 Jun 2023
Update on 08 Jun 2023
@author: Hanyu Jin
version: 0.1
"""
from numpy import ndarray
from numba import njit
@njit()
def find_compound_event(event_times: list[list[tuple]], max_gap: ndarray = None):
    """
    Calculate the conditional probabilities and the cascading event periods for the given event time ranges.

    :param event_times: A list of lists, where each sublist contains tuples representing the start and end times of
    events for a variable. For example: [(start, end), ...]
    :param max_gap: The maximum allowed time gap (in time
    points) between consecutive events for them to be considered as a chain (cascade). :return: A list containing the
    conditional probability and the cascade periods for each pair of adjacent events.
    """
    result = []
    m = len(event_times)  # Number of variables

    # Iterate over each pair of adjacent event lists (e.g., ENSO -> SPEI, SPEI -> WILDFIRE)
    for i in range(m - 1):
        cause_events = event_times[i]  # Current variable's event time periods
        effect_events = event_times[i + 1]  # Next variable's event time periods

        cause_count = 0  # Count of events in the cause
        cause_and_effect_count = 0  # Count of events where both cause and effect happen
        cascade_periods = []  # Store the periods where cascading events occur

        # Check for overlap or time gap between the events of the current variable and the next variable
        for cause_start, cause_end in cause_events:
            for effect_start, effect_end in effect_events:
                # Check if the start of the effect event is within 2 time points before or after the end of the cause
                # event
                if (cause_end - max_gap[i] <= effect_start <= cause_end + max_gap[i]) or \
                        (cause_start - max_gap[i] <= effect_end <= cause_start + max_gap[i]):  # Cascade happens
                    cause_count += 1
                    cause_and_effect_count += 1
                    # Record the time period where cascading events happen
                    overlap_start = max(cause_start, effect_start)
                    overlap_end = min(cause_end, effect_end)
                    cascade_periods.append((overlap_start, overlap_end))

        # calculate conditional probability
        if cause_count == 0:  # If no cause events happened, the probability is 0
            conditional_probability = 0.0
        else:
            conditional_probability = cause_and_effect_count / cause_count

        # append the conditional probability and cascade periods to the result
        result.append([conditional_probability, cascade_periods])

    return result


@njit()
def find_consecutive(seq: list, delta: int, max_gap_length: int, max_gap_count: int) -> list[tuple]:
    """
    Find consecutive sequences of True values in a binary sequence, considering gaps between the sequences.

    :param seq: A list of boolean values (True/False) representing the sequence to be analyzed.
    :param delta: Minimum length of the consecutive True sequence to be considered valid.
    :param max_gap_length: Maximum allowed length of consecutive False values (gaps) within a sequence.
    :param max_gap_count: Maximum number of consecutive False values allowed (the number of breaks) within a sequence.
    :return: A list of tuples, where each tuple contains the start and end index of a valid consecutive sequence.
    """
    result = []  # List to store the valid consecutive sequences as tuples of (start, end)
    n = len(seq)  # Length of the input sequence
    i = 0  # Index pointer to traverse through the sequence

    # Traverse through the sequence
    while i < n:
        if seq[i]:  # If the current element is True (1)
            start = i  # Mark the start of the current subsequence
            broken = 0  # Counter to keep track of the number of breaks (False values)
            gap_length = 0  # Counter for the length of the current gap of False values
            last_true_index = i  # The index of the last encountered True value

            while i < n:
                if seq[i]:  # If the current element is True
                    last_true_index = i  # Update the last True index
                    gap_length = 0  # Reset the gap length as we're encountering a True value
                else:  # If the current element is False
                    gap_length += 1  # Increment the gap length
                    if gap_length > max_gap_length:  # If the gap length exceeds the allowed maximum
                        break  # Break the loop, as this sequence is no longer valid
                    broken += 1  # Increment the number of breaks (False values encountered)
                    if broken > max_gap_count:  # If the number of breaks exceeds the allowed maximum
                        break  # Break the loop, as this sequence is no longer valid

                # If the gap length or the number of breaks exceeds the limits, break the inner loop
                if gap_length > max_gap_length or broken > max_gap_count:
                    break

                i += 1  # Move to the next index

            end = last_true_index  # The index of the last True value in the current sequence
            result.append((start, end))  # Add the current valid subsequence to the result as a tuple (start, end)
        else:
            i += 1  # Skip the False value and move to the next index

    # Filter the results to only include subsequences that meet the minimum length requirement (delta)
    result_filter = [x for x in result if (x[1] - x[0] + 1) >= delta]
    return result_filter  # Return the filtered list of consecutive subsequences




def apply_operator(a, b, op):
    """
    this function convert original signal to bool array
    :param a: original signal
    :param b: threshold
    :param op: operator of lower and upper bounds
    :return: binary array
    """
    operations = {
        ('ge', 'le'): lambda a, b: (a >= b[0]) & (a <= b[1]),
        ('g', 'le'): lambda a, b: (a > b[0]) & (a <= b[1]),
        ('ge', 'l'): lambda a, b: (a >= b[0]) & (a < b[1]),
        ('g', 'l'): lambda a, b: (a > b[0]) & (a < b[1])
    }

    if op in operations:
        return operations[(op)](a, b)
    else:
        raise ValueError(f"Unsupported operation combination: {op}, the valid is ('ge', 'le'), ('g', 'le'), "
                         f"('ge', 'l'), ('g', 'l')")

# add two significant test possion and surrogate methods
