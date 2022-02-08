#!/usr/bin/env ruby
module IOUtils
  #
  # num_comma(n)
  #   e.g. n = 10000 
  #    num_comma(n) => '10,000'
  #
  def num_comma(n)
    nn = ""
    n.to_s.reverse.split("").each_with_index do |v, i|
      nn += "," if i % 3 == 0 && i > 0
      nn += v
    end
    return nn.reverse
  end

  def print_blue(str,type='b')
    ct_print_str(str,'blue',type)
  end

  def print_red(str,type='b')
    ct_print_str(str,'red',type)
  end

  def print_green(str,type='b')
    ct_print_str(str,'green',type)
  end

  def print_cyan(str,type='b')
    ct_print_str(str,'cyan',type)
  end

  def print_black(str,type='b')
    ct_print_str(str,'black',type)
  end

  def print_bold(str)
    print "#{ct_escape "1"}#{str}#{ct_reset}"
  end

  def print_underline(str)
    print "#{ct_escape "4"} #{str} #{ct_reset}"
  end

  def ct_print_str(str,color='blue',type='b')
    case(type.strip.downcase)
    when /^u$/,/^b$/ then
      print "#{ct_color(color,type.strip.downcase)}#{str}#{ct_reset}"
    else
      print "#{str}"
    end
  end
  def ct_color(color,type)
    ansi_color = {
      'black'  => 30,
      'red'    => 31,
      'green'  => 32,
      'yello'  => 33,
      'blue'   => 34,
      'magenda' => 35,
      'cyan'   => 36,
      'white'   => 37
    }
    begin
      unless ansi_color[color.downcase].nil? then
        num = ansi_color[color.downcase]
      else
        num = 30
      end
      case(type.downcase)
      when /^u$/
        ct_underline num
      when /^b$/
        ct_bold num
      end
    rescue
      print_red "error was happen.\n"
      exit 1
    end
  end

  def ct_reset
    ct_escape(0)
  end

  def ct_bold(n)
    ct_escape("1;#{n}")
  end

  def ct_underline(n)
    ct_escape("4;#{n}")
  end

  def ct_escape(n)
    "\033[#{n}m" if $stdout.tty?
  end
end
