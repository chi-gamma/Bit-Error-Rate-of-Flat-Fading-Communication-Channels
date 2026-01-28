import random
from collections import deque
from dataclasses import dataclass
import heapq



@dataclass(frozen=True)
class RollingResult:
    mean: float
    variance: float
    minimum: float
    maximum: float



class RollingStatistics:
    def __init__(self, window_size, ddof=1):
        if window_size <= ddof:
            raise ValueError("Window size must be greater than degrees of freedom")
        self.window_size = window_size
        self.ddof = ddof
        self.window = deque()
        self.min_window = deque()
        self.max_window = deque()
        self.sum_x = 0.0
        self.sum_x2 = 0.0
        self.index = 0
        self.result = {'mean' : float('nan'), 'variance' : float('nan'), 'minimum' : float('nan')}

    def update(self, price):
        # minimum logic (window is monotonic increasing)
        while self.min_window and self.min_window[-1][0] >= price:
            self.min_window.pop()

        self.min_window.append((price, self.index))
        
        while self.min_window[0][1] < self.index - self.window_size + 1:
            self.min_window.popleft()
            
            
        # maximum logic (window is monotonic decreasing)
        while self.max_window and self.max_window[-1][0] <= price:
            self.max_window.pop()

        self.max_window.append((price, self.index))
        
        while self.max_window[0][1] < self.index - self.window_size + 1:
            self.max_window.popleft()
        
        self.index += 1
        
        # mean and variance logic     
        self.window.append(price)
        self.sum_x += price
        self.sum_x2 += price * price
        
        if len(self.window) > self.window_size:
            old = self.window.popleft()
            self.sum_x -= old
            self.sum_x2 -= old * old

        return self._calculate()
            
    def _calculate(self):
        n = self.window_size
        if len(self.window) < n:
            result = RollingResult(
                mean = float('nan'),
                variance = float('nan'),
                minimum = float('nan'),
                maximum = float('nan')
                )            
            return result # not enough data
        
        mean = self.sum_x / n
        ssd = self.sum_x2 - ( self.sum_x * self.sum_x ) / n # sum of square deviations
        var = ssd / (n - self.ddof) if ssd > 0.0 else 0.0
        
        result = RollingResult(
            mean = mean,
            variance = var,
            minimum = self.min_window[0][0],
            maximum = self.max_window[0][0]
            )

        return result
        

rs = RollingStatistics(30)
N = 1000
prices = [ random.uniform(a=5020.0, b=6200.0) for _ in range(1000) ]

means = []
variances = []
minimums = []
maximums = []

for price in prices:
    res = rs.update(price)
    m, v, low, high = res.mean, res.variance, res.minimum, res.maximum
    means.append(m)
    variances.append(v)
    minimums.append(low)
    maximums.append(high)
    






# check
import numpy as np 
prices = np.array(prices)
means_ch = []
variances_ch = []
minimums_ch = []
maximums_ch = []
for i in range(30, len(prices)+1):
    arr = prices[i-30 : i]
    means_ch.append(arr.mean())
    variances_ch.append(arr.var(ddof=1))
    minimums_ch.append(arr.min())
    maximums_ch.append(arr.max())
    




    


# import heapq
# from collections import defaultdict, deque
# import math

# class RollingMedian:
#     def __init__(self, window_size: int):
#         self.w = window_size
#         self.low = []   # max-heap via negatives
#         self.high = []  # min-heap
#         self.delayed = defaultdict(int)
#         self.low_size = 0
#         self.high_size = 0
#         self.window = deque()  # store values to know what expires

#     def _prune(self, heap, is_low: bool):
#         # Remove heap tops that are marked for deletion
#         while heap:
#             x = -heap[0] if is_low else heap[0]
#             if self.delayed[x] > 0:
#                 heapq.heappop(heap)
#                 self.delayed[x] -= 1
#                 if self.delayed[x] == 0:
#                     del self.delayed[x]
#             else:
#                 break

#     def _rebalance(self):
#         # Ensure low has either same number of valid elems as high, or one more
#         if self.low_size > self.high_size + 1:
#             self._prune(self.low, True)
#             x = -heapq.heappop(self.low)
#             self.low_size -= 1
#             heapq.heappush(self.high, x)
#             self.high_size += 1
#         elif self.low_size < self.high_size:
#             self._prune(self.high, False)
#             x = heapq.heappop(self.high)
#             self.high_size -= 1
#             heapq.heappush(self.low, -x)
#             self.low_size += 1

#         self._prune(self.low, True)
#         self._prune(self.high, False)

#     def update(self, x: float):
#         # Add new value
#         self.window.append(x)

#         if not self.low:
#             heapq.heappush(self.low, -x)
#             self.low_size += 1
#         else:
#             self._prune(self.low, True)
#             if x <= -self.low[0]:
#                 heapq.heappush(self.low, -x)
#                 self.low_size += 1
#             else:
#                 heapq.heappush(self.high, x)
#                 self.high_size += 1

#         # Expire old value if window too large
#         if len(self.window) > self.w:
#             out = self.window.popleft()
#             self.delayed[out] += 1

#             # Decide which heap it belongs to using current boundary
#             self._prune(self.low, True)
#             if out <= -self.low[0]:
#                 self.low_size -= 1
#                 if out == -self.low[0]:
#                     self._prune(self.low, True)
#             else:
#                 self.high_size -= 1
#                 self._prune(self.high, False)

#         self._rebalance()

#         # Return median if we have full window
#         if len(self.window) < self.w:
#             return math.nan

#         if self.w % 2 == 1:
#             self._prune(self.low, True)
#             return -self.low[0]
#         else:
#             self._prune(self.low, True)
#             self._prune(self.high, False)
#             return (-self.low[0] + self.high[0]) / 2.0
    







class RollingMedian:
    def __init__(self, window_size):
        if window_size <= 0:
            raise ValueError("window size must be positive")

        self.delayed = dict()
        self.low = []
        self.high = []
        self.window = deque()
        self.low_size = 0
        self.high_size = 0
        self.window_size = window_size
        self.median_fn = self._median_even if window_size % 2 == 0 else self._median_odd()
        
        
    def _median_even(self):
        return (self.low_max() + self.high_min()) / 2.0
    
    def _median_odd(self):
        return self.low_max()
    
    def low_max(self) -> float:
        """ Returns the current maximum value in the low heap"""
        return -self.low[0]
    
    def high_min(self) -> float:
        """ Returns the current minimum value in the high heap"""
        return self.high[0]
    
    
    def clean_low(self) -> None:
        """ Removes expired elements from the top of the low heap """
        while self.low:
            x = -self.low[0]
            if self.delayed.get(x, 0) > 0:
                heapq.heappop(self.low)
                self.delayed[x] -= 1
                if self.delayed[x] == 0:
                    del self.delayed[x]
            else:
                break
        
    def clean_high(self) -> None:
        """ Removes expired elements from the top of the high heap """
        while self.high:
            x = self.high[0]
            if self.delayed.get(x, 0) > 0:
                heapq.heappop(self.high)
                self.delayed[x] -= 1
                if self.delayed[x] == 0:
                    del self.delayed[x]
            else:
                break
        
    def balance(self) -> None:
        """ Maintain size equilibrium (low_size == high_size OR low_size == high_size + 1) """
        while self.low_size > self.high_size + 1:
            self.clean_low() # ensure you are not transferring an expired price from low heap
            x = -heapq.heappop(self.low)
            self.low_size -= 1
            heapq.heappush(self.high, x)
            self.high_size += 1
            self.clean_low() # ensure the new maximum of low heap is not an expired price
            
        while self.low_size < self.high_size:
            self.clean_high() # ensure you are not transferring an expired price from high heap
            x = heapq.heappop(self.high)
            self.high_size -= 1
            heapq.heappush(self.low, -x)
            self.low_size += 1
            self.clean_high() # ensure the new maximum of high heap is not an expired price
            
        
    def update(self, x: float) -> None:
        self.window.append(x)
        
        self.clean_low()
        if not self.low or x <= self.low_max():
            heapq.heappush(self.low, -x)
            self.low_size += 1
        else:
            heapq.heappush(self.high, x)
            self.high_size += 1
            
        if len(self.window) > self.window_size:
            old = self.window.popleft()
            # flag this old price as expired
            self.delayed[old] += 1
            
            self.clean_low()
            self.clean_high()
            if self.low and old <= self.low_max():
                self.low_size -= 1
            else:
                self.high_size -= 1
                
        self.balance()
        
        if len(self.window) < self.window_size:
            return float('nan')
        
        self.clean_low()
        self.clean_high()
        return self.median_fn()
        
        
        
rm = RollingMedian(30)
N = 1000
prices = [ random.uniform(a=5020.0, b=6200.0) for _ in range(1000) ]
medians = []
for price in prices:
    res = rm.update(price)
    medians.append(res)
    
    
    
    
    m, v, low, high = res.mean, res.variance, res.minimum, res.maximum
    means.append(m)
    variances.append(v)
    minimums.append(low)
    maximums.append(high)
        
        
        
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    



    


