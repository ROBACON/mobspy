from mobspy import *

if __name__ == '__main__':

    Age, Color, Mortal = BaseSpecies(3)
    Age.young, Age.old
    Color.red, Color.green, Color.blue

    Mortal >> Zero [1]

    A,B,C = New(Age*Color*Mortal, 3)

    with Any.old.red :
        A >> B [1]
        B >> C [1]
        C >> A [1]

    with Any.young.green :
        A(10)
        B(10)
        C(10)

    S = Simulation(A | B |C)

    S.method = 'deterministic'
    S.duration = 15
    S.output_event = True

    with Any.old, S.event_condition(A <= 1) :
        with Color.red :
            A(20)
            C(5)

        with Any.green :
            A(10)
            C(10)

    S.run()