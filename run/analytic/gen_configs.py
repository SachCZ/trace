import functools
import itertools
import yaml


def generate_configs(grad_types, segments, factors):
    base_config = {
        "mesh": {
            "x0": 0,
            "x1": 0.75,
            "y0": 0,
            "y1": 1.1107207345395915,
            "x_segments": 8,
            "y_segments": 8,
            "randomize_factor": 0,
        },
        "ele_dens_profile": "1.7145e+24/2 * (1 + x^2)",
        "rays_filename": "",
        "lasers": [{
            "wavelength": 25.5e-7,
            "direction": {"x": 1, "y": 1, },
            "start_point": {"x": -0.0001, "y": 0.6},
            "end_point": {"x": -0.0001, "y": -0.0001},
            "rays_count": 50,
            "grad_type": ""
        }]
    }

    for grad_type, segment_count, factor in itertools.product(grad_types, segments, factors):
        config = base_config.copy()
        config["lasers"][0]["grad_type"] = grad_type
        if grad_type == "mfem":
            config["lasers"][0]["grad_boundary_x"] = "1.7145e+24 * x"
            config["lasers"][0]["grad_boundary_y"] = "0"

        if segment_count == 32 and factor == 0:
            config["mesh_output_filename"] = "./output/mesh{}{}{}.mesh".format(segment_count, grad_type, factor)

        config["mesh"]["x_segments"] = segment_count
        config["mesh"]["y_segments"] = segment_count
        config["mesh"]["randomize_factor"] = factor
        config["rays_filename"] = "./output/rays{}{}{}.msgpack".format(segment_count, grad_type, factor)

        filename = "input/config{}{}{}.yaml".format(segment_count, grad_type, factor)
        with open(filename, "w") as f:
            yaml.dump(config, f)

        yield filename


def main():
    grad_types = ["ls", "integral", "mfem"]
    segments = [8, 16, 32, 64, 128, 256, 512]
    factors = [0, 0.02, 0.04]
    config_names = list(generate_configs(grad_types, segments, factors))
    print(functools.reduce(lambda res, name: res + "\n" + name, config_names))


if __name__ == '__main__':
    main()
