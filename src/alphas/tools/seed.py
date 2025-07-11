import random as rand

if __name__ == "__main__":
    seed_length = 9
    seed = ""
    for _ in range(seed_length): seed += str(rand.randint(0,9))
    print("seed:",seed)
