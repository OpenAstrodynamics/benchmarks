import icatt.elements
import icatt.kepler
import icatt.lambert
import icatt.dopri

times = 100000

if __name__ == "__main__":
    icatt.elements.benchmark(times)
    icatt.kepler.benchmark(times)
    icatt.lambert.benchmark(times)
    icatt.dopri.benchmark(times)
