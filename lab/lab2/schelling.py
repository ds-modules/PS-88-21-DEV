# import matplotlib
# matplotlib.use("TkAgg")
from enum import Enum
from random import choices, sample
from collections import namedtuple, deque
import matplotlib.pyplot as plt
from itertools import product, chain

Pref = namedtuple('Pref', ['in_frac', 'color'])

GroupWeights = namedtuple(
    'GroupWeights', ['friendly', 'shy', 'fearful', 'racist'])


def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(
        zip(handles, labels)) if l not in labels[:i]]
    plt.legend(*zip(*unique), loc='upper right')


class Group(Enum):
    # Tuple of desired ingroup fraction, plotting color
    FRIENDLY = Pref(0, 'b')  # No need for any change
    SHY = Pref(.25, 'g')  # Want at least a quarter
    FEARFUL = Pref(.5, 'r')  # Want half
    RACIST = Pref(.8, 'y')  # Nearly only ingroup

    @staticmethod
    def random_pick(weights: GroupWeights = None):
        """
        Pick a random group type
        """
        return choices(list(Group), weights)[0]


class Agent:
    def __str__(self):
        return f'{self.x}, {self.y}, {self.group}'

    def __init__(self, group: Group, x: int, y: int):
        self.x, self.y, self.group = x, y, group
        # self.satified = None

    def plot(self, ax):
        plt.scatter(self.x, self.y, c=self.group.value.color,
                    label=self.group.name)


class Grid:

    def __init__(self, N: int, vacancy_rate: float, weights: GroupWeights = None):
        self.N = N
        self.moved_any = True
        self.matrix = [[None]*N for _ in range(N)]
        self.fig, self.ax = plt.subplots(figsize=(10, 10))
        for i, j in product(range(N), repeat=2):
            self.matrix[i][j] = Agent(Group.random_pick(weights), i, j)
        for i, j in product(sample(range(N), int(N*vacancy_rate)), repeat=2):
            self.matrix[i][j] = None

    def vacancies(self) -> [(int, int)]:
        return [(x, y) for x, y in product(range(self.N), repeat=2) if self.matrix[x][y] is None]

    def neighbors(self, a):
        """
        Return the neighbors of the agent at the location
        """
        return [self.matrix[n_x][n_y] for n_x, n_y in product(*(range(n-1, n+2) for n in (a.x, a.y))) if (n_x, n_y) != (a.x, a.y) and all(0 <= n < self.N for n in (n_x, n_y)) and self.matrix[n_x][n_y] is not None]

    def satisfied(self, a):
        """
        If agent a has a sufficient number of the ingroup surrounding them, they are satisfied
        """
        neighbors = self.neighbors(a)
        return len(list(filter(lambda n: n.group == a.group, neighbors))) >= len(neighbors)*a.group.value.in_frac

    def plot(self):
        for i, j in product(range(self.N), repeat=2):
            if self.matrix[i][j] is not None:
                self.matrix[i][j].plot(self.ax)
        legend_without_duplicate_labels(self.ax)
        # self.fig.show()
        plt.show()

    def move(self, a: Agent):
        if a is not None and not self.satisfied(a):
            acceptable_vacancies = [
                (x, y) for x, y in self.vacancies() if self.satisfied(Agent(a.group, x, y))]
            if len(acceptable_vacancies) > 0:
                spot_x, spot_y = acceptable_vacancies[0]
                mid_x, mid_y = (spot_x+a.x)/2, (spot_y+a.y)/2
                plt.annotate("moving", (a.x, a.y),
                             xytext=(mid_x, mid_y), arrowprops={'arrowstyle': '->'})
                plt.annotate("moving", (spot_x, spot_y),
                             xytext=(mid_x, mid_y), arrowprops={'arrowstyle': '->'})

                self.matrix[a.x][a.y] = None
                a.x, a.y = spot_x, spot_y
                self.matrix[spot_x][spot_y] = a
                return True
        return False

    def move_first(self):
        agents = deque(filter(lambda a: a is not None,
                              chain.from_iterable(self.matrix)))
        moved_any = None
        while (moved_any is None or not moved_any) and len(agents) > 0:
            a = agents.popleft()
            moved_any = self.move(a)
        return moved_any

    def run(self):
        moved_any = None
        while moved_any is None or moved_any:
            # print(moved_any)
            self.plot()
            moved_any = self.move_first()
        self.plot()


if __name__ == "__main__":
    m = Grid(10, .4, GroupWeights(.1, .3, 0.4, .2))
    # print(m)
    # m.plot()
    m.run()
    # m.schelling()
    m.plot()
    # plt.show()
