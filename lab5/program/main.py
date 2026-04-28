import cmath
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VARIANT_FILE = ROOT / "task" / "var.md"
EPS = 1e-10


def parse_variant(path):
    text = path.read_text(encoding="utf-8")

    variant_match = re.search(r"variant:\s*(\d+)", text)
    values_match = re.search(r"values:\s*(.+)", text)
    if not variant_match or not values_match:
        raise ValueError("File var.md must contain 'variant:' and 'values:' lines.")

    variant = int(variant_match.group(1))
    values = [float(item.replace(",", ".")) for item in values_match.group(1).split()]
    return variant, values


def dft(values):
    n = len(values)
    result = []
    for k in range(n):
        total = 0j
        for j, value in enumerate(values):
            total += value * cmath.exp(-2j * math.pi * k * j / n)
        result.append(total)
    return result


def inverse_dft(amplitudes):
    n = len(amplitudes)
    result = []
    for j in range(n):
        total = 0j
        for k, value in enumerate(amplitudes):
            total += value * cmath.exp(2j * math.pi * k * j / n)
        result.append(total / n)
    return result


def fft_decimation_by_time(values):
    n = len(values)
    if n % 2 != 0:
        raise ValueError("The number of samples must be even for radix-2 decimation.")

    even_values = values[0::2]
    odd_values = values[1::2]
    even_amplitudes = dft(even_values)
    odd_amplitudes = dft(odd_values)

    amplitudes = [0j] * n
    for k in range(n // 2):
        multiplier = cmath.exp(-2j * math.pi * k / n)
        amplitudes[k] = even_amplitudes[k] + multiplier * odd_amplitudes[k]
        amplitudes[k + n // 2] = even_amplitudes[k] - multiplier * odd_amplitudes[k]

    return even_values, odd_values, even_amplitudes, odd_amplitudes, amplitudes


def inverse_fft_decimation_by_time(even_amplitudes, odd_amplitudes):
    even_values = inverse_dft(even_amplitudes)
    odd_values = inverse_dft(odd_amplitudes)

    values = []
    for even_value, odd_value in zip(even_values, odd_values):
        values.append(even_value)
        values.append(odd_value)

    return even_values, odd_values, values


def clean_number(value, digits=6):
    if abs(value) < EPS:
        value = 0.0
    return f"{value:.{digits}f}"


def format_complex(value, digits=6):
    real = 0.0 if abs(value.real) < EPS else value.real
    imag = 0.0 if abs(value.imag) < EPS else value.imag

    if imag == 0.0:
        return clean_number(real, digits)
    if real == 0.0:
        return f"{clean_number(imag, digits)}i"

    sign = "+" if imag >= 0 else "-"
    return f"{clean_number(real, digits)} {sign} {clean_number(abs(imag), digits)}i"


def print_real_vector(title, values):
    print(title)
    print(f"{'j':>3} {'value':>18}")
    for j, value in enumerate(values):
        number = value.real if isinstance(value, complex) else value
        print(f"{j:>3} {clean_number(number):>18}")
    print()


def print_complex_vector(title, values, index_name="k"):
    print(title)
    print(f"{index_name:>3} {'value':>28}")
    for index, value in enumerate(values):
        print(f"{index:>3} {format_complex(value):>28}")
    print()


def print_verification(title, original, restored, index_name="j", original_name="original", restored_name="restored"):
    print(title)
    print(f"{index_name:>3} {original_name:>14} {restored_name:>28} {'abs error':>14}")
    max_error = 0.0
    for j, (source, value) in enumerate(zip(original, restored)):
        error = abs(source - value)
        max_error = max(max_error, error)
        print(
            f"{j:>3} "
            f"{format_complex(source):>14} "
            f"{format_complex(value):>28} "
            f"{error:>14.3e}"
        )
    print(f"max error = {max_error:.3e}")
    print()


def main():
    variant, values = parse_variant(VARIANT_FILE)

    print("=" * 78)
    print("LAB 5: DISCRETE FOURIER TRANSFORM")
    print("=" * 78)
    print(f"Variant: {variant}")
    print(f"Number of samples: n = {len(values)}")
    print()

    print_real_vector("Initial signal y(j):", values)

    print("-" * 78)
    print("TASK 1. Direct DFT and inverse DFT")
    print("-" * 78)
    amplitudes = dft(values)
    restored = inverse_dft(amplitudes)

    print_complex_vector("DFT amplitudes F(k):", amplitudes)
    print_complex_vector("Inverse DFT result y(j):", restored, index_name="j")
    print_verification("DFT reconstruction check:", values, restored)

    print("-" * 78)
    print("TASK 2. Fast Fourier transform by time decimation")
    print("-" * 78)
    (
        even_values,
        odd_values,
        even_amplitudes,
        odd_amplitudes,
        fft_amplitudes,
    ) = fft_decimation_by_time(values)

    print_real_vector("Even-index samples y^0(j) = y(2j):", even_values)
    print_real_vector("Odd-index samples y^1(j) = y(2j + 1):", odd_values)
    print_complex_vector("DFT by even samples F^0(k):", even_amplitudes)
    print_complex_vector("DFT by odd samples F^1(k):", odd_amplitudes)
    print_complex_vector("FFT amplitudes F(k):", fft_amplitudes)

    even_restored, odd_restored, fft_restored = inverse_fft_decimation_by_time(
        even_amplitudes,
        odd_amplitudes,
    )

    print_complex_vector("Inverse transform by even samples y^0(j):", even_restored, index_name="j")
    print_complex_vector("Inverse transform by odd samples y^1(j):", odd_restored, index_name="j")
    print_complex_vector("Restored signal from inverse FFT blocks y(j):", fft_restored, index_name="j")

    print_verification(
        "FFT amplitudes compared with direct DFT:",
        amplitudes,
        fft_amplitudes,
        index_name="k",
        original_name="direct DFT",
        restored_name="FFT",
    )
    print_verification("FFT reconstruction check:", values, fft_restored)

    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()
