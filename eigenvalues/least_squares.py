from random import randint
from random import random
from random import uniform

max_iter = 1000

def fit(xy):
	x = xy[0]
	y = xy[1]
	return (2*x + 3*y - 5)**2 + (x + 3*y - 4)**2 + (2*x + 3*y - 4)**2

def generate_combinations(xy_list, N):
	length = len(xy_list) - 1
	for i in range(N):
		rand_x = randint(0, length)
		rand_y = randint(0, length)
		yield [xy_list[rand_x][0], xy_list[rand_y][1]]

def get_random():
	return uniform(0,2)

def generate_random_xy(N):
	for i in range(N):
		x = get_random()
		y = get_random()
		yield [x, y]

def mutate(xy, mutation_prob):
	if random() < mutation_prob:
		xy[0] = get_random()
	if random() < mutation_prob:
		xy[1] = get_random()
	return xy

def genetic(sample_size=100, mutation_prob=0.1, keep_best=30):
	sample = list(generate_random_xy(100))
	offspring_size = sample_size - keep_best

	mutation_prob_ = mutation_prob
	prev_fitness = 0

	k = 0
	while True:
		k += 1
		fitness = [(xy, fit(xy)) for xy in sample]
		fitness.sort(key=lambda x: x[1], reverse=False)

		if k % 10 == 0:
			print("{}: Best pair ({}, {}) with fitness {}" \
				.format(k, fitness[0][0][0], fitness[0][0][1], fitness[0][1]))
			if fitness[0][1] == prev_fitness:
				mutation_prob_ = min(1.2 * mutation_prob_, 0.9)
			else:
				mutation_prob_ = mutation_prob
			prev_fitness = fitness[0][1]

		sample = [x[0] for x in fitness[:keep_best]]

		for x in generate_combinations(sample, offspring_size):
			sample.append(mutate(x, mutation_prob_))

if __name__ == "__main__":
	genetic(mutation_prob=0.5, keep_best=400, sample_size=20000)

