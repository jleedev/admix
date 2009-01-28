class AccountsController < ApplicationController
  def login
    @title = "Log in"
    if request.post?
      @current_user = User.find_by_name params[:login]
      if @current_user
        session[:user_id] = @current_user.id
        redirect_to :controller => :marker_file
      else
        flash[:error] = "Incorrect username"
      end
    end
  end
end
