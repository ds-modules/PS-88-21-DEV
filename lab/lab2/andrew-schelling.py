from math import sqrt
from copy import deepcopy
from operator import itemgetter
from enum import Enum
import matplotlib.pyplot as plt
from collections import namedtuple
Pref = namedtuple('Pref', ['in_weight', 'out_weight'])


class AgentType(Enum):
    A = 'A'
    B = 'B'
    EMPTY = '_'

    def __str__(self):
        return self.value

    @staticmethod
    def other(agent):
        assert agent.agent_type != AgentType.EMPTY
        return AgentType.A if agent.agent_type == AgentType.B else AgentType.B


class Agent:
    @staticmethod
    def swap(agent_x, agent_y):
        mid_x = (agent_x.position[0]+agent_y.position[0])/2
        plt.annotate("moving", (*agent_x.position, 0), xytext=(mid_x,
                                                               -0.004), arrowprops={'arrowstyle': '->'})
        plt.annotate("moving", (*agent_y.position, 0), xytext=(mid_x,
                                                               -0.004), arrowprops={'arrowstyle': '->'})
        temp = deepcopy(agent_x.position)
        agent_x.position = agent_y.position
        agent_y.position = temp

    def __init__(self, coords: int, agent_type: AgentType):
        self.position = coords
        self.agent_type = agent_type

    @classmethod
    def from_other(a, coords, other):
        return a(coords, other.agent_type)

    def __str__(self):
        return f"Agent of type {self.agent_type} at position {self.position} with preferences {self.preferences}"


class Map:
    @staticmethod
    def distance(agent1, agent2):
        assert len(agent1.position) == len(agent2.position)
        return sqrt(sum(((agent1.position[i] - agent2.position[i])**2 for i in range(len(agent1.position)))))

    @staticmethod
    def from_string(agents_string):
        agents = []
        for idx, agent_type in enumerate(agents_string):
            agents.append(Agent([idx], AgentType(agent_type)))
        return Map(agents)

    def __init__(self, agents: list):
        self.agents = agents

    def __str__(self):
        map_str = "Map with agents: \n"
        for a in self.agents:
            map_str += "Agent at position " + str(a.position) + "\n"
        return map_str

    def plot(self):
        # print(self.get_coords())
        for a in self.agents:
            color = 'r' if a.agent_type == AgentType.A else 'g' if a.agent_type == AgentType.B else 'b'
            plt.scatter(a.position, 0, c=color)
            plt.annotate(str(a.agent_type), (*a.position, 0), xytext=(*a.position,
                                                                      0.001), rotation=60)
        plt.show()

    def get_coords(self):
        return [a.position for a in self.agents]

    def get_type(self, agent_type):
        return list(filter(lambda a: a.agent_type == agent_type, self.agents))

    def dist_to_type(self, agent, agent_type):
        others = self.get_type(agent_type)
        others = filter(lambda other_agent: other_agent is not agent, others)
        return min((Map.distance(agent, other_agent) for other_agent in others))
    # def avoid_other_type(self, agent, other_type):
    #     # Swap it with empty locales further from the other types
    #     swappable = self.get_type(AgentType.EMPTY)
    #     furthest_swappable = max(
    #         swappable, key=lambda swap_agent: self.dist_to_type(swap_agent, other_type))
    #     if self.dist_to_type(agent, other_type) < self.dist_to_type(furthest_swappable, other_type):
    #         # print('hey', agent.agent_type)
    #         print(
    #             f"swapping agents at {agent.position}, {furthest_swappable.position}""")
    #         Agent.swap(agent, furthest_swappable)
    #         return True
    #     return False

    def pos_score(self, agent):
        # Score a position based on its closeness to its preferred type and distance from other type (higher is better)
        return -1*agent.pref.in_weight*self.dist_to_type(agent, agent.agent_type) + agent.pref.out_weight*self.dist_to_type(agent, AgentType.other(agent))

    def a_score(self, idx):
        current_score = 0
        current_score -= 1 if idx > 0 and self.agents[idx-1] == AgentType.B
        current_score -= 1 if idx < len(self.agents) - 1 and self.agents[idx+1] == AgentType.B
        return current_score

    def b_score(self, idx):
        current_score = 0
        current_score += 1 if idx > 0 and self.agents[idx-1] != AgentType.EMPTY
        current_score += 1 if idx < len(self.agents) - 1 and self.agents[idx+1] != AgentType.EMPTY
        return current_score

    def vacancy_score(self, agent: Agent, vacancy_pos):
        return self.pos_score(Agent.from_other(vacancy_pos, agent))

    def move(self):
        for idx, agent in enumerate(self.agents):
            vacant = self.get_type(AgentType.EMPTY)
            if agent.agent_type != AgentType.EMPTY:
                if agent.agent_type == 'A':
                    vacant_scores = max(vacant, key=lambda v: self.a_score(agent, v.position))
        # vacant_weights = {s: self.pos_score(
    def min_other(self, a: Agent):

    def move(self, agent: Agent):
        metric = self.min_other if agent.agent_type == AgentType.A else self.max_all

        vacant = self.get_type(AgentType.EMPTY)
        best_vacant = max(
            vacant, key=lambda v: self.vacancy_score(agent, v.position))
        # vacant_weights = {s: self.pos_score(
        #     s, agent.agent_type, other_type) for s in vacant}
        # best_swappable = max(vacant_weights.items(),
        #                      key=itemgetter(1))[0]
        # best_open = max(vacant, key=lambda a: self.pos_score(
        #     Agent(a.position, agent.agent_type, agent.weights), agent.agent_type, other_type))
        # print(best_swappable == self.pos_score(best_open, ))
        # print(self.weighted_preference(agent, agent.agent_type, other_type), agent.prefences)
        # print(vacant_weights)
        if self.pos_score(agent) < self.vacancy_score(agent, best_vacant.position):
            print(self.pos_score(agent), self.vacancy_score(
                agent, best_vacant.position))
            print(
                f"swapping agents at {agent.position}, {best_vacant.position}""")
            Agent.swap(agent, best_vacant)
            return True
        return False

    def move_all(self):
        agents = sorted(self.agents, key=lambda a: a.position)
        moved_any = False
        for a in filter(lambda a: a.agent_type != AgentType.EMPTY, agents):
            a_moved = self.move(a)
            if a_moved:
                moved_any = True
                self.plot()
                # return moved_any
        return moved_any

    def schelling(self):
        runtime = 0
        moved_any = True
        self.plot()
        while moved_any:
            some_moved = self.move_all()
            if not some_moved:
                moved_any = False
                plt.title('Finished')
                self.plot()
            runtime += 1
            if runtime > 100:
                print("There may not be a stable arrangement")
                break


if __name__ == "__main__":
    m = Map.from_string('ABABAB___', Pref(0, 1), Pref(1, -1))
    # print(m.agents[0].preferences)
    # print(m)
    # m.plot()
    m.schelling()
    # m.plot()
