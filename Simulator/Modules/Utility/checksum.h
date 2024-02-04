#pragma once

int checksum(int* data, size_t length) {
	int sum = 0;
	for (size_t i = 0; i < length; i++) {
		sum += data[i];
	}
	return sum;
}