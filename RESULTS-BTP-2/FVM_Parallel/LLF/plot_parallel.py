import matplotlib.pyplot as plt

def read_grid(file: str) -> list:
    """
        Function to read the finite volume computational grid
    """
    with open(file, "r") as f:
        grid = list(map(float, f.readline().split(" ")[:-1]))

    return grid

def read_density(file: str) -> list[list]:
    """
        Function to read the density 
        -> first element of the returned list is initial density
        -> last element of the returned list is final density
    """
    with open(file, "r") as f:
        lines = f.readlines()

    density = []

    for line in lines:
        density.append( list(map(float, line.split(" ")[:-1])) )

    return density

def plot_density(grid: list, density: list[list], **kwargs) -> None:
    """
        Function to plot the initial and final density 
    """
    # plt.plot(grid, density[0], label="Initial Density", color="black", linestyle="--")
    plt.plot(grid, density[-1], label="Final Density", color="blue")

    plt.title(kwargs["MAIN_TITLE"])
    plt.legend()

    plt.savefig(f"PLOTS/PARALLEL-{kwargs['MAIN_TITLE']}.png")

    # plt.show()

def main(MODE: str) -> None:
    """
        main function
    """
    # PROBLEM = "Moving Contact Wave"
    # PROBLEM = "Blast Wave Problem"
    # PROBLEM = "LAX Problem"
    PROBLEM = "TORO-1"
    # PROBLEM = "TORO-2"
    # PROBLEM = "TORO-3"

    FILE_GRID = f"parallel_grid.txt"
    FILE_DENSITY = f"parallel_density.txt"

    # read the computational grid and the density
    grid = read_grid(FILE_GRID)
    density = read_density(FILE_DENSITY)

    # Plot the density initially and finally
    plot_density(grid, density, MAIN_TITLE=PROBLEM+"(PARALLEL)")
    
if __name__ == "__main__":
    main("NORMAL")
    # main("REF")




