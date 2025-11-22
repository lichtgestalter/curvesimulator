
framesss =            [9, 10, 11,  11,  11,  11,  11,  11]
total_iterationsss = [10, 10, 10, 100, 109, 110, 111, 112]

for frames, total_iterations in zip(framesss, total_iterationsss):
    print(f"{total_iterations=}   {frames=}   {(total_iterations - 1) // frames + 1=}")