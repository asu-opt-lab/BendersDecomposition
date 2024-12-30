using BendersDecomposition

I = 50
J = 50
S = 256

r1 = 3
r2 = 5
r3 = 10

for i in 1:5
    data = generate_stochastic_capacited_facility_location(I,J,S,r1)
    write_stochastic_capacited_facility_location_problem(data; filename="f$(I)-c$(J)-s$(S)-r$(r1)-$(i).json")
end

for i in 1:5
    data = generate_stochastic_capacited_facility_location(I,J,S,r2)
    write_stochastic_capacited_facility_location_problem(data; filename="f$(I)-c$(J)-s$(S)-r$(r2)-$(i).json")
end

for i in 1:5
    data = generate_stochastic_capacited_facility_location(I,J,S,r3)
    write_stochastic_capacited_facility_location_problem(data; filename="f$(I)-c$(J)-s$(S)-r$(r3)-$(i).json")
end