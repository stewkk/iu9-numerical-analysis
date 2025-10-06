{ config, pkgs, ... }: {
  environment.systemPackages = [pkgs.gnuplot pkgs.julia-bin pkgs.unstable.cursor-cli];
  # ...
}
