module AdmixWeb

  def self.prettify str
    str.split("\n").map do |line|
      "<tr><td>" + line.gsub(/\ +/,"<td>")
    end.join("\n")
  end

end
